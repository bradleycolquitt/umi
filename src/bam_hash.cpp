#include <bam_hash.h>

using namespace std;

void UmiHash::update(BamRecord* record)
{
    pair<unordered_map<int, int>::const_iterator, bool> result;
    uint32_t umi = record->get_umi();
    result = umi_map.emplace(umi, 1);
    if (!result.second) {
        ++umi_map[umi];
    }
}
void BarcodeHash::update(BamRecord * record)
{
    pair<unordered_map<int, unique_ptr<UmiHash> >::const_iterator, bool> result;
    result = umi_map.emplace(record->get_bc(), unique_ptr<UmiHash>(new UmiHash));
    result.first->second->update(record);
}

void PositionHash::update(BamRecord * record)
{
    pair<unordered_map<int, unique_ptr<BarcodeHash> >::const_iterator, bool> result;
    result = barcode_map.emplace(record->get_pos(), unique_ptr<BarcodeHash>(new BarcodeHash));
    result.first->second->update(record);
}

// Main Bamhash constructor
BamHash::BamHash(const char* bam_fname, const char* anno_fname, const char* dest_fname, const char* barcodes_fname, int umi_length, int bc_min_qual)
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
    //ofstream out(dest_path.string(), ofstream::out);
    outfile.open(dest_path.string(), ofstream::out);
    // tmp_path = dest_path;
    // tmp_path += ".tmp";

    // merge_path = dest_path.parent_path();
    // merge_path /= "merge.db";



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

BamHash::~BamHash() {
    // sqlite3_close(conns[0]);
    // sqlite3_close(conns[1]);
    sam_close(bam);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    for (auto& position : position_map)
    {
        //delete position.second;
    }

    //remove_tmp_files();
}

// Read in barcodes file and load sequences into barcodes vector
void BamHash::set_barcodes(const char* fname, vector<vector<int> >& vec_p) {
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

void BamHash::insert_anno(const string & read_id, const string & gene_id) {
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

void BamHash::update_maps(BamRecord * record)
{

    pair<unordered_map<string, unique_ptr<PositionHash>  >::const_iterator, bool> result;
    result = position_map.emplace(record->get_gene_id(), unique_ptr<PositionHash>(new PositionHash()));
    result.first->second->update(record);
}

void hash_annotation(BamHash* bamhash)
{
    ifstream anno(bamhash->get_anno_fname());
    string read_id, assigned, gene_id;
    const string assigned_val("Assigned");
    while(anno >> read_id >> assigned >> gene_id) {
        if (assigned_val.compare(assigned) == 0) {
            bamhash->insert_anno(read_id, gene_id);
        }
    }
}

void process_read1(BamHash* bamhash, BamRecord * record, bam1_t* b)
{
    int used_offset = 0;
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

void BamHash::print_results() {
    //for (auto& position : bamhash->get_position_map())
    for (auto& gene : position_map)
    {
        unordered_map<int, unique_ptr<BarcodeHash> >::iterator position = gene.second->begin();
        for (; position != gene.second->end(); ++position)
        //for (auto& barcode : position.second )
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
    hts_itr_t* bam_itr;
    for (int tid = 0; tid < header->n_targets; ++tid)
    {
        BamRecord * record = new BamRecord(tid);
        bam_itr = bam_itr_queryi(idx, tid, 0, get_rlen(tid));
        hash_reads_tid(this, record, bam_itr);
        delete record;
    }
    bam_itr_destroy(bam_itr);
}

// void BamHash::run() {
//     clock_t start;
//     start = std::clock();

//     //try {
//     hash_reads();
//     print_time(start);

//     // }  catch (runtime_error &e) {
//     //     cerr << " ! Runtime error: " << e.what() << endl;
//     //     exit(1);
//     // }
//}
