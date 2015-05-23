#include <bam_filter.h>

/* Main Bamhash constructor */
BamHash::BamHash(const char* bam_fname, const char* anno_fname, const char* dest_fname,                  const char* barcodes_fname, int umi_length, int bc_min_qual)
    : bam_fname(bam_fname)
    , anno_fname(anno_fname)
    , dest_fname(dest_fname)
    , bc_min_qual(bc_min_qual)
    {
    bam = sam_open(bam_fname, "rb");
    header = sam_hdr_read(bam);
    idx = bam_index_load(bam_fname);

    char ** chrom_names = header->target_name;
    for (int i = 0 ; i < header->n_targets ; ++i) {
        chroms[i] = chrom_names[i];
    }

    boost::filesystem::path dest_fname_path(dest_fname);
    dest_path = boost::filesystem::complete(dest_fname_path);
    outbam = sam_open(dest_path.string().c_str(), "wb");
    bam_hdr_write(outbam->fp.bgzf, header);
    //outfile.open(dest_path.string(), ofstream::out);

    //Barcode setup
    char bc_path[255];
    strcpy(bc_path, BC_PATH);
    strcat(bc_path, barcodes_fname);
    strcat(bc_path, ".txt");
    try {
        set_barcodes(bc_path, barcodes);
    } catch (exception &e) {
        cout << e.what() << endl;
    }

    //Positions of UMI and barcodes relative to 5' end of read
    sequence_pos.push_back(0);
    sequence_pos.push_back(umi_length-1);
    sequence_pos.push_back(umi_length);
    size_t bc_length = barcodes[0].size();
    sequence_pos.push_back(sequence_pos[2] + bc_length - 1);

    // defines offsets used during barcode search
    int offsets_array[] = {0,-1,1};
    bc_offsets.assign(offsets_array, offsets_array + sizeof(offsets_array) / sizeof(int));
}

BamHash::~BamHash()
{
    sam_close(bam);
    sam_close(outbam);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
}

// Read in barcodes file and load sequences into barcodes vector
void BamHash::set_barcodes(const char* fname, vector<vector<int> >& vec_p)
{
    ifstream bc_s(fname, ifstream::in);
    if (!bc_s.good()) {
         throw runtime_error("Error: Barcode file not found.");
    }

    string line;
    while(getline(bc_s, line)) {
        if ((bc_s.rdstate() & ifstream::failbit ) != 0) {
           cout << "Error" << endl;
        }
        vector<string> sline = split(line, '\t');
        vector<int> s2i = seq2int(sline[1]);
        vec_p.push_back(s2i);
    }
}

bool UmiHash::update(BamRecord* record)
{
    pair<unordered_map<int, int>::const_iterator, bool> result;
    uint32_t umi = record->get_umi();
    result = umi_map.emplace(umi, 1);
    if (!result.second) {
        ++umi_map[umi];
    }
    return result.second;
}

bool BarcodeHash::update(BamRecord * record)
{
    pair<unordered_map<int, unique_ptr<UmiHash> >::const_iterator, bool> result;
    result = umi_map.emplace(record->get_bc(), unique_ptr<UmiHash>(new UmiHash));
    return result.first->second->update(record);
}

bool PositionHash::update(BamRecord * record)
{
    pair<unordered_map<int, unique_ptr<BarcodeHash> >::const_iterator, bool> result;
    result = barcode_map.emplace(record->get_pos(), unique_ptr<BarcodeHash>(new BarcodeHash));
    return result.first->second->update(record);
}

// bool StrandHash::update(BamRecord * record)
// {
//     pair<unordered_map<int, unique_ptr<PositionHash> >::const_iterator, bool> result;
//     result = position_map.emplace(record->get_pos(), unique_ptr<PositionHash>(new PositionHash));
//     return result.first->second->update(record);
// }

bool BamHash::update_maps(BamRecord * record)
{
    pair<unordered_map<bool, unique_ptr<PositionHash> >::const_iterator, bool> result;
    result = strand_position_map.emplace(record->get_strand(),
                                  unique_ptr<PositionHash>(new PositionHash()));
    return result.first->second->update(record);
}

void print_bam_record(BamHash * bamhash, bam1_t* b)
{
    bam_write1(bamhash->get_outbam()->fp.bgzf, b);
}

void process_read1(BamHash* bamhash, BamRecord * record, bam1_t* b)
{
    int used_offset = 0;
    if (bad_cigar(b)) return;
    if (filter_multi_reads(b)) return;

    record->set_strand(b);
    record->set_position(b->core.pos);
    record->set_bc(bamhash, b, &used_offset);
    record->set_umi(bamhash, b, used_offset);
    bool found = bamhash->update_maps(record);
    if (! found )
    {
        print_bam_record(bamhash, b);
    }
}

int hash_reads_tid(BamHash* bamhash, BamRecord * record, hts_itr_t* bam_itr)
{
    bam1_t* b = bam_init1();
    int result;
    while ((result = sam_itr_next(bamhash->get_bam(), bam_itr, b)) >= 0)
    {
        //if (!bamhash->find_read(b, record)) continue;
        if ((b->core.flag&BAM_FREAD1) != 0)
        {
            process_read1(bamhash, record, b);
        }
    }
    bam_destroy1(b);
    return 0;
}

void BamHash::hash_reads()
{
    cout << "Hashing reads..." << endl;
    hts_itr_t* bam_itr;

    int total = header->n_targets;
    int measure = round(total)/10;
    for (int tid = 0; tid < header->n_targets; ++tid)
    {
        if (measure>1)
        {
            if ((tid % measure) == 0)
            {
                cout << "-";
                cout.flush();
            }
        }
        BamRecord * record = new BamRecord(tid);
        bam_itr = bam_itr_queryi(idx, tid, 0, get_rlen(tid));
        hash_reads_tid(this, record, bam_itr);
        delete record;
        clear_map();
        bam_itr_destroy(bam_itr);
    }
    cout << endl;

}

void BamHash::clear_map()
{
    strand_position_map.clear();
}
