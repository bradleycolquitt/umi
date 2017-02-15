#include <bam_hash2.h>

using namespace std;

/**************************
    Update map functions
***************************/
bool UmiHash::update(shared_ptr<BamRecord> record)
{
    //bool print = false;
    pair<unordered_map<string, int>::iterator, bool> result;
    string umi = record->get_umi2();

    result = umi_map.emplace(umi, 1);

    if (!result.second) {
        /* UMI already present in map so increment */
        //if (print) cout << "initial: " << result.first->second << endl;
        ++(result.first->second);
        //if (print) cout << "update: " << result.first->second << endl;
    }
    return result.second;
}

bool BarcodeHash::update(shared_ptr<BamRecord> record)
{
    pair<unordered_map<int, unique_ptr<UmiHash> >::const_iterator, bool> result;
    result = umi_map.emplace(record->get_bc(), unique_ptr<UmiHash>(new UmiHash));
    return result.first->second->update(record);
}

bool PositionHash::update(shared_ptr<BamRecord> record)
{
    pair<unordered_map<int, unique_ptr<BarcodeHash> >::const_iterator, bool> result;
    result = barcode_map.emplace(record->get_pos(), unique_ptr<BarcodeHash>(new BarcodeHash));
    return result.first->second->update(record);
}

/****************************
   BamHash functions
*****************************/
BamHash::BamHash(const char* bam_fname, const char* fastq_fname, const char* anno_fname,   \
                 const char* dest_fname, const char* barcodes_fname,    \
                 int umi_length,                                        \
                 int bc_min_qual,                                       \
                 int i5, int i7, bool to_txt)
    : bam_fname(bam_fname)
    , fastq_fname(fastq_fname)
    , anno_fname(anno_fname)
    , dest_fname(dest_fname)
    , bc_min_qual(bc_min_qual)
    , i5(i5)
    , i7(i7)
    , to_txt(to_txt)
    {

    /* Check existance of files*/
    // if (!std::ifstream(fastq_fname))
    // {
    //     std::cout << "! Fastq file (" << fastq_fname << ") doesn't exist!" << std::endl;
    //     exit(1);
    // }

    // if (!std::ifstream(anno_fname))
    // {
    //     std::cout << "! featureCounts file (" << anno_fname << ") doesn't exist!" << std::endl;
    //     exit(1);
    // }

    // boost::filesystem::path dest_fname_path(dest_fname);
    // dest_path = boost::filesystem::complete(dest_fname_path);

    if (to_txt)
    {
        //outfile.open(dest_path.string(), ofstream::out);
    }
    else // setup database
    {
        try
        {
            const char* dest_path_c = dest_path.string().c_str();
            open_connection(dest_path_c, &conn, true);
            sqlite3_busy_timeout(conn, 1000);
        }
        catch (sql_exception &e)
        {
            cerr << "! Error: opening databases, " << e.what() << endl;
            exit(1);
        }

        try
        {
            exec_multithread(conn, "CREATE TABLE IF NOT EXISTS counts (i5 int, i7 int, bc int, gene_id text, umi int, count int);" );
            //execute(conn, "CREATE TABLE IF NOT EXISTS counts (i5 int, i7 int, bc int, gene_id text, umi int, count int);");
        }
        catch (sql_exception &e)
        {
            cerr << "! Error: creating table, " << e.what() << endl;
        }
    }

    //Barcode setup
    char bc_path[255];
    strcpy(bc_path, BC_PATH);
    strcat(bc_path, barcodes_fname);
    strcat(bc_path, ".txt");
    try {
        set_barcodesA(bc_path, barcodesA);
    } catch (exception &e) {
        cout << e.what() << endl;
    }

    //Positions of UMI and barcodes relative to 5' end of read
    sequence_pos.push_back(0);
    sequence_pos.push_back(umi_length-1);
    sequence_pos.push_back(umi_length);
    size_t bc_length = sizeof(barcodes[0]) / sizeof(barcodes[0][0]);
    sequence_pos.push_back(sequence_pos[2] + bc_length - 1);

    // defines offsets used during barcode search
    int offsets_array[] = {0,-1,1};
    bc_offsets.assign(offsets_array, offsets_array + sizeof(offsets_array) / sizeof(int));
}

BamHash::~BamHash()
{
   int result = sqlite3_close(conn);
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

// As set_barcodes but read in as vector of strings
void BamHash::set_barcodesA(const char* fname, vector<const char*>& vec_p)
{
    ifstream bc_s(fname, ifstream::in);
    cout << fname << endl;
    if (!bc_s.good()) {
         throw runtime_error("Error: Barcode file not found.");
    }

    string line;
    while(getline(bc_s, line)) {
        if ((bc_s.rdstate() & ifstream::failbit ) != 0) {
           cout << "Error" << endl;
        }
        std::istringstream iss(line);
        string id;
        string bc;
        if (iss >> id >> bc)
        {
            const char* bc_s = bc.c_str();
            //vec_p.push_back(strdup(bc.c_str()));
            vec_p.push_back(bc_s);
        }
    }
}

void BamHash::insert_anno(const string & read_id, const string & gene_id)
{
    qname_bamrecord.emplace(read_id, shared_ptr<BamRecord>(new BamRecord(read_id, gene_id)));
}

void BamHash::insert_bam_pos(string & qname, bam1_t* br)
{
    bam_map.emplace(qname, shared_ptr<BamRecord>(new BamRecord(br->core.tid, br->core.pos)));
}

bool BamHash::find_read(bam1_t * b, shared_ptr<BamRecord> & record)
{
    string read_id = string(bam_get_qname(b));
    unordered_map<string, shared_ptr<BamRecord> >::iterator result;
    if ((result = qname_bamrecord.find(read_id)) != qname_bamrecord.end())
    {
        record = result->second;
        return true;
    }
    else
    {
        return false;
    }
}

bool BamHash::find_read(char * qname, shared_ptr<BamRecord> & record)
{
    string read_id = string(qname);
    unordered_map<string, shared_ptr<BamRecord> >::iterator result;
    if ((result = qname_bamrecord.find(read_id)) != qname_bamrecord.end())
    {
        record = result->second;
        return true;
    }
    else
    {
        return false;
    }
}

bool BamHash::update_maps(shared_ptr<BamRecord> record)
{
    pair<unordered_map<string, unique_ptr<BarcodeHash> >::const_iterator, bool> result;
    result = position_map.emplace(record->get_gene_id(),
                                  unique_ptr<BarcodeHash>(new BarcodeHash()));
    result.first->second->update(record);
    return true;
}

void hash_annotation(BamHash* bamhash)
{
    cout << endl;
    cout << "----Loading annotation----" << endl;
    ifstream anno(bamhash->get_anno_fname());

    string line;
    string read_id, assigned, gene_id;
    const string assigned_val("Assigned");
    int hash_count = 0;
    while (getline(anno, line))
    {
        istringstream iss(line);
        if (!(iss >> read_id >> assigned >> gene_id)) continue;
        if (assigned_val.compare(assigned) == 0) {
            hash_count++;
            bamhash->insert_anno(read_id, gene_id);
        }
    }
    cout << "Number of reads in annotation set: " << hash_count << endl;
    cout << endl;
}

void hash_bam(BamHash* bamhash)
{
    cout << "here" << endl;

    samFile* in_bam = sam_open(bamhash->get_bam(), "r");
    bamhash->set_bam(in_bam);

    bam_hdr_t* header = sam_hdr_read(bamhash->get_bam_file());

    string bam_fname1 = string(bamhash->get_bam());
    string bam_idx_name = bam_fname1 + ".bai";

    if (!boost::filesystem::exists(bam_idx_name))
    {
        cout << "Building BAM index..." << endl;
        bam_index_build(bamhash->get_bam(), 0);
    }

    hts_idx_t* idx = bam_index_load(bam_idx_name.c_str());
    bamhash->set_bam_idx(idx);
    bam1_t* br = bam_init1();

    // loop through tids
    // iterate through reads and add to unordered_map
    int tid = 0;
    int res = 0;
    const char* qname;
    string qname_s;

    cout << "here" << endl;
    while (tid < header->n_targets)
    {
        cout << tid << endl;

        hts_itr_t* itr = bam_itr_queryi(idx, tid, 0, header->target_len[tid]);
        //hts_itr_t* itr = bam_itr_queryi(idx, tid, 0, 50);


        //cout << "here" << endl;
        //res = bam_itr_next(in_bam, itr, br);
        //cout << res << endl;
        while((res = bam_itr_next(in_bam, itr, br)) >= 0) {
             qname = bam_get_qname(br);
             qname_s = string(qname);
             bamhash->insert_bam_pos(qname_s, br);
             //bam_map->emplace(qname, shared_ptr<bam1_t*>(&br));
         }
        tid++;
    }

}
void parse_fastq(BamHash* bamhash) {
    gzFile fp;
    kseq_t* seq;
    int l;

    fp = gzopen(bamhash->get_fastq(), "r");
    seq = kseq_init(fp);

    int used_offset = 0;
    shared_ptr<BamRecord> record;
    //int complete = 0;
    while ((l = kseq_read(seq)) >= 0)
    {
        if (!bamhash->find_read(seq->name.s, record)) continue; // sets record to qname_bamrecord
        record->set_bc(bamhash, seq->seq.s, seq->qual.s, &used_offset);
        record->set_umi(bamhash, seq->seq.s, seq->qual.s, used_offset);
        used_offset = 0;
        bamhash->update_maps(record);
    }
    kseq_destroy(seq);
    gzclose(fp);
}

void parse_bam_by_fastq(BamHash* bamhash) {
    gzFile fp;
    kseq_t* seq;
    int l;

    fp = gzopen(bamhash->get_fastq(), "r");
    seq = kseq_init(fp);

    int used_offset = 0;
    shared_ptr<BamRecord> record;
    //int complete = 0;
    while ((l = kseq_read(seq)) >= 0)
    {
        if (!bamhash->find_read(seq->name.s, record)) continue;
        record->set_bc(bamhash, seq->seq.s, seq->qual.s, &used_offset);
        bamhash->split_bam_by_record(record);
        used_offset = 0;

    }
    kseq_destroy(seq);
    gzclose(fp);
}

void BamHash::split_bam_by_record(shared_ptr<BamRecord> record)
{

    hts_itr_t* itr = bam_itr_queryi(bam_idx, record->get_tid(), record->get_pos(), record->get_pos() + 1);
    bam1_t* br = bam_init1();
    int res = 0;
    const char* qname;
    while ((res = bam_itr_next(in_bam, itr, br))>=0)
    {
        qname = bam_get_qname(br);
        const char * read_id_c = record->get_read_id().c_str();
        if (strcmp(qname, read_id_c))
        {
            unordered_map<int,samFile*>::const_iterator got = out_bam_map.find (record->get_bc());

            if ( got == out_bam_map.end() )
            {
                continue;
            }
            else
            {
                std::cout << got->first << " is " << got->second;
                bam_write1(got->second->fp.bgzf, br);
            }
        }
    }

}

void BamHash::print_results()
{
    vector<string> file_header;
    file_header.push_back("i5");
    file_header.push_back("i7");
    file_header.push_back("bc");
    file_header.push_back("gene_id");
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
        unordered_map<int, unique_ptr<UmiHash> >::iterator barcode = gene.second->begin();
        for (; barcode != gene.second->end(); ++barcode)
        {
            unordered_map<string, int>::iterator umi = barcode->second->begin();
            for (; umi != barcode->second->end(); ++umi)
            {
                outfile << this->i5 << "\t"
                     << this->i7 << "\t"
                     << barcode->first << "\t"
                     << gene.first << "\t"
                     << umi->first << "\t"
                     << umi->second << endl;
            }
        }
    }
}

void BamHash::write_to_db()
{
    const char* tail = 0;
    char insert_sql[BUFFER_SIZE];
    sqlite3_stmt * insert_stmt;
    sprintf(insert_sql, "INSERT INTO counts VALUES (@I5, @I7, @BC, @GEN, @UMI, @COUNT);");
    sqlite3_prepare_v2(conn, insert_sql, BUFFER_SIZE, \
                       &insert_stmt, &tail);

    execute(conn, "BEGIN");
    for (auto& gene : position_map)
    {
        unordered_map<int, unique_ptr<UmiHash> >::iterator barcode = gene.second->begin();
        for (; barcode != gene.second->end(); ++barcode)
        {
            unordered_map<string, int>::iterator umi = barcode->second->begin();
            for (; umi != barcode->second->end(); ++umi)
            {
                sqlite3_bind_int(insert_stmt, 1, this->i5);
                sqlite3_bind_int(insert_stmt, 2, this->i7);
                sqlite3_bind_int(insert_stmt, 3, barcode->first);
                sqlite3_bind_text(insert_stmt, 4, gene.first.c_str(), -1, SQLITE_TRANSIENT);
                //sqlite3_bind_int(insert_stmt, 5, position->first);
                //sqlite3_bind_text(insert_stmt, 6, umi->first.c_str(), -1, SQLITE_TRANSIENT);
                //sqlite3_bind_int(insert_stmt, 7, umi->second);
                sqlite3_bind_text(insert_stmt, 5, umi->first.c_str(), -1, SQLITE_TRANSIENT); /* UMI */
                sqlite3_bind_int(insert_stmt, 6, umi->second); /* count */

                int result = 0;
                //sqlite3_step(insert_stmt);
                //if ((result = sqlite3_step(insert_stmt)) == SQLITE_BUSY) {
                if ((result = sqlite3_step(insert_stmt)) == SQLITE_BUSY )
                {
                    throw sql_exception(result, sqlite3_errmsg(conn));
                }
                sqlite3_clear_bindings(insert_stmt);
                sqlite3_reset(insert_stmt);
            }
        }
    }
    execute(conn, "COMMIT");
    sqlite3_finalize(insert_stmt);
}

void BamHash::split_bam()
{
    hash_bam(this);
    setup_out_bams();
    parse_bam_by_fastq(this);
}

void BamHash::setup_out_bams()
{
    /* set up BAMs for all libraries */
    string bam_fname1 = string(bam_fname);
    samFile* in_bam = sam_open(bam_fname, "r");
    size_t lastindex = bam_fname1.find_last_of(".");
    string rawname = bam_fname1.substr(0, lastindex);
    string basename = boost::filesystem::basename(bam_fname1);
    string outdir = rawname;

    boost::filesystem::path dir(outdir);

    if(!(boost::filesystem::exists(dir))){
        if (boost::filesystem::create_directory(dir))
            std::cout << "....Successfully Created !" << std::endl;
    }

    /* create unordered_map of bam_objects keyed by bc */
    //unordered_map<string, samFile* > bam_map;
    bam_hdr_t* header = sam_hdr_read(in_bam);
    int bc_index = 1;
    for (auto it = barcodesA.begin(); it != barcodesA.end(); ++it)
    {
       const char* bc_cur = *it;
       string fname_cur = outdir + "/" + basename + "_bc" + to_string(bc_index) + ".bam";
       const char* fname_cur_cs = convert_to_cstr(fname_cur);
       samFile* out_bam = sam_open(fname_cur_cs, "wb");
       cout << bc_cur << fname_cur << endl;
       out_bam_map.insert({bc_index, out_bam});
       bc_index++;
    }

    // write header
    for (auto it = out_bam_map.begin(); it != out_bam_map.end(); ++it)
    {
       cout << it->first << endl;
       bam_hdr_write(it->second->fp.bgzf, header);
       sam_close(it->second);
    }
}

void BamHash::write_bam()
{
    //seq =  kseq_init(f_rd);
    // bam = sam_open(bam_fname, "r");

    /* set up BAMs for all libraries */
    string bam_fname1 = string(bam_fname);
    samFile* in_bam = sam_open(bam_fname, "r");
    size_t lastindex = bam_fname1.find_last_of(".");
    string rawname = bam_fname1.substr(0, lastindex);
    string basename = boost::filesystem::basename(bam_fname1);
    string outdir = rawname;

    boost::filesystem::path dir(outdir);

    if(!(boost::filesystem::exists(dir))){
        if (boost::filesystem::create_directory(dir))
            std::cout << "....Successfully Created !" << std::endl;
    }

    /* create unordered_map of bam_objects keyed by bc */
    unordered_map<string, samFile* > bam_map;
    bam_hdr_t* header = sam_hdr_read(in_bam);
    for (auto it = barcodesA.begin(); it != barcodesA.end(); ++it)
    {
       const char* bc_cur = *it;
       string fname_cur = outdir + "/" + basename + "_" + bc_cur + ".bam";
       const char* fname_cur_cs = convert_to_cstr(fname_cur);
       samFile* out_bam = sam_open(fname_cur_cs, "wb");
       cout << bc_cur << fname_cur << endl;
       bam_map.insert({bc_cur, out_bam});
    }

    // write header
    for (auto it = bam_map.begin(); it != bam_map.end(); ++it)
    {
       cout << it->first << endl;
       bam_hdr_write(it->second->fp.bgzf, header);
       sam_close(it->second);
    }

    //gzFile fp;
    //kseq_t* seq;
    //int l;

    //fp = gzopen(fastq_fname, "r");
    //seq = kseq_init(fp);

    //int used_offset = 0;
    //shared_ptr<BamRecord> record;
    //int complete = 0;

    //const hts_idx_t* idx = bam_index_load(bam_fname);
    //hts_itr
    //bam1_t* br = bam_init1();

    //const char* fastq_name;
    //const char* qname;
    // iterate through fastq and bam
    // while ((l = kseq_read(seq)) >= 0 && (res = bam_itr_next(in_bam, itr, br)) >= 0)
    // {
    //     fastq_name = seq->name.s;
    //     qname = bam_get_qname(br);
    //     cout << fastq_name << " " << qname << endl;
    //     //hts_itr_t* read_itr = bam_itr_querys(idx, header, seq->name.s);

    //     //int res = bam_readrec(in_bam->fp.bgzf, NULL, read, read_itr->tid, read_itr->beg, read_itr->end);
    //     //cout << read->core.tid << endl;
    //     //if (!bamhash->find_read(seq->name.s, record)) continue; // sets record to qname_bamrecord
    //     //record->set_bc(bamhash, seq->seq.s, seq->qual.s, &used_offset);
    //     //record->set_umi(bamhash, seq->seq.s, seq->qual.s, used_offset);
    //     //used_offset = 0;
    //     //bamhash->update_maps(record);
    // }
    // kseq_destroy(seq);
    // gzclose(fp);


    // get barcode
    // get bam1
    // write bam read to approriate file

    // for ( auto it = qname_bamrecord.begin(); it != qname_bamrecord.end(); ++it )
    // {
    // //     bam1_t *q = bam_init1();
    // }

    //     printf("%s, %s\n", seq->seq.s, seq->qual.s);

    //     //`q->data` structure: qname-cigar-seq-qual-aux

    //     q->l_data = seq->name.l+1+(int)(1.5*seq->seq.l+(seq->seq.l % 2 != 0)); // +1 includes the tailing '\0'
    //     if (q->m_data < q->l_data) {
    //         q->m_data = q->l_data;
    //         kroundup32(q->m_data);
    //         q->data = (uint8_t*)realloc(q->data, q->m_data);
    //     }
    //     q->core.flag = BAM_FMUNMAP;
    //     q->core.l_qname = seq->name.l+1; // +1 includes the tailing '\0'
    //     q->core.l_qseq = seq->seq.l;
    //     q->core.n_cigar = 0; // we have no cigar sequence
    //     memcpy(q->data, seq->name.s, q->core.l_qname); // first set qname
    //     uint8_t *s = bam_get_seq(q);
    //     for (int i = 0; i < q->core.l_qseq; ++i){

    //         bam1_seq_seti(s, i, seq_nt16_table[seq->seq.s[i]]);
    //     }
    //     s = bam_get_qual(q);
    //     for (int i = 0; i < seq->qual.l; ++i) {
    //         s[i] = seq->qual.s[i];
    //     }
    //     dump_read(q);
    //     int ret = bam_write1(fp, q);
    //     printf("return value: %d\n", ret);
    //     bam_destroy1(q);
    // }
    // kseq_destroy(seq); // free seq
    // gzclose(f_rd); // close fastq file
    // bgzf_close(fp); // close bam file


    // // validate bam file
    // BGZF *fp_new = bgzf_open(fn_out,"r");
    // bam1_t *q = bam_init1();
    // int ret = bam_read1(fp_new, q);
    // printf("%d\n", ret);
    // printf("%d, %d\n", (char *)bam_get_aux(q) - (char *)q->data , q->l_data);
    // dump_read(q);
    // return 0;
}
void BamHash::hash_reads()
{
    cout << "----Parsing information reads----" << endl;
    parse_fastq(this);
    cout << endl;

}
