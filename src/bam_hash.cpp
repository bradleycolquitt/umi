#include <bam_hash.h>

using namespace std;

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

// Main Bamhash constructor
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
    outfile.open(dest_path.string(), ofstream::out);

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

void BamHash::insert_anno(const string & read_id, const string & gene_id)
{
    anno_map.emplace(read_id, gene_id);
}

bool BamHash::find_read(bam1_t * b, BamRecord * record)
{
    record->set_read_id(bam_get_qname(b));
    unordered_map<string,string>::const_iterator result;

    if ((result = anno_map.find(record->get_read_id())) != anno_map.end())
    {
        record->set_gene_id(result->second);
        return true;
    } else {
        return false;
    }
}

bool BamHash::update_maps(BamRecord * record)
//void BamHash::update_maps(BamRecord * record)
{
    pair<unordered_map<string, unique_ptr<PositionHash> >::const_iterator, bool> result;
    result = position_map.emplace(record->get_gene_id(),
                                  unique_ptr<PositionHash>(new PositionHash()));
    result.first->second->update(record);
    return true;
}

void hash_annotation(BamHash* bamhash)
{
    cout << "Loading annotation..." << endl;
    ifstream anno(bamhash->get_anno_fname());
    string line;
    string read_id, assigned, gene_id;
    const string assigned_val("Assigned");
//    const string test_val("augustus_masked-scaffold_0-processed-gene-4.2");
    //while(anno >> read_id >> assigned >> gene_id) {
    while (getline(anno, line))
    {
        istringstream iss(line);
        if (!(iss >> read_id >> assigned >> gene_id)) continue;
        if (assigned_val.compare(assigned) == 0) {
//            if (test_val.compare(gene_id) == 0) cout << read_id << "\t" << gene_id << endl;
            bamhash->insert_anno(read_id, gene_id);
        }
    }
}

void process_read1(BamHash* bamhash, BamRecord * record, bam1_t* b)
{
    int used_offset = 0;
    if (bad_cigar(b)) return;
    if (filter_multi_reads(b)) return;

    record->set_position(b->core.pos);
    record->set_bc(bamhash, b, &used_offset);
    record->set_umi(bamhash, b, used_offset);
    bamhash->update_maps(record);
}

int hash_reads_tid(BamHash* bamhash, BamRecord * record, hts_itr_t* bam_itr)
{
    bam1_t* b = bam_init1();
    int result;
    while ((result = sam_itr_next(bamhash->get_bam(), bam_itr, b)) >= 0) {
        if (!bamhash->find_read(b, record)) continue;
        if ((b->core.flag&BAM_FREAD1) != 0) {
            process_read1(bamhash, record, b);
        }
    }
    bam_destroy1(b);
    return 0;
}

void BamHash::print_results()
{
    vector<string> file_header;
    file_header.push_back("bc");
    file_header.push_back("gene_id");
    file_header.push_back("position");
    file_header.push_back("umi");
    file_header.push_back("count");

    vector<string>::iterator iter = file_header.begin();
    outfile << *iter;
    ++iter;
    for (; iter != file_header.end() ; ++iter)
    {
        outfile << "\t" << *iter;
    }
    outfile << endl;

    for (auto& gene : position_map)
    {
        unordered_map<int, unique_ptr<BarcodeHash> >::iterator position = gene.second->begin();
        for (; position != gene.second->end(); ++position)
        {
            unordered_map<int, unique_ptr<UmiHash> >::iterator barcode = position->second->begin();
            for (; barcode != position->second->end(); ++barcode)
            {
                unordered_map<int, int>::iterator umi = barcode->second->begin();
                for (; umi != barcode->second->end(); ++umi)
                {
                    outfile << barcode->first << "\t"
                         << gene.first << "\t"
                         << position->first << "\t"
                         << umi->first << "\t"
                         << umi->second << endl;
                }
            }
        }
    }
}

void BamHash::hash_reads()
{
    cout << "Hashing reads..." << endl;
    hts_itr_t* bam_itr;

    int total = header->n_targets;
    int measure = round(total)/10;
    for (int tid = 0; tid < header->n_targets; ++tid)
    {
        //cout << tid << endl;
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
        bam_itr_destroy(bam_itr);
    }
    cout << endl;

}
