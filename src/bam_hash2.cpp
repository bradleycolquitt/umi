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
BamHash::BamHash(const char* fastq_fname, const char* anno_fname,       \
                 const char* dest_fname, const char* barcodes_fname,    \
                 int umi_length,                                        \
                 int bc_min_qual,                                       \
                 int i5, int i7, bool to_txt)
    : fastq_fname(fastq_fname)
    , anno_fname(anno_fname)
    , dest_fname(dest_fname)
    , bc_min_qual(bc_min_qual)
    , i5(i5)
    , i7(i7)
    , to_txt(to_txt)
    {

    /* Check existance of files*/
    if (!std::ifstream(fastq_fname))
    {
        std::cout << "! Fastq file (" << fastq_fname << ") doesn't exist!" << std::endl;
        exit(1);
    }

    if (!std::ifstream(anno_fname))
    {
        std::cout << "! featureCounts file (" << anno_fname << ") doesn't exist!" << std::endl;
        exit(1);
    }

    boost::filesystem::path dest_fname_path(dest_fname);
    dest_path = boost::filesystem::complete(dest_fname_path);

    if (to_txt)
    {
        outfile.open(dest_path.string(), ofstream::out);
    }
    else // setup database
    {
        try
        {
            open_connection(dest_path.string().c_str(), &conn, true);
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
            vec_p.push_back(strdup(bc.c_str()));
        }
    }
}

void BamHash::insert_anno(const string & read_id, const string & gene_id)
{
    qname_bamrecord.emplace(read_id, shared_ptr<BamRecord>(new BamRecord(read_id, gene_id)));
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

void parse_fastq(BamHash* bamhash) {
    gzFile fp;
    kseq_t* seq;
    int l;

    fp = gzopen(bamhash->get_fastq(), "r");
    seq = kseq_init(fp);

    int used_offset = 0;
    shared_ptr<BamRecord> record;
    int complete = 0;
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

void BamHash::hash_reads()
{
    cout << "----Parsing information reads----" << endl;
    parse_fastq(this);
    cout << endl;

}
