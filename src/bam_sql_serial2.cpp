#include <bam_sql_serial2.h>
//#include <bam_utils.h>
//#include <string_utils.h>
//#include <sqlite3.h>
//#include <boost/regex.hpp>
//#include <sstream>

using namespace std;

#ifdef DEBUG
#define DEBUG_LOG(x) cerr << x << endl;
#else
#define DEBUG_LOG(x)
#endif

// Main BamDB constructor
BamDB::BamDB(const char* bam_fname, const char* dest_fname, const char* barcodes_fname, int umi_length, int bc_min_qual, bool paired_end)
    : bam_fname(bam_fname)
    , dest_fname(dest_fname)
    , paired_end(paired_end)
    , bc_min_qual(bc_min_qual)
    {
        bam = sam_open(bam_fname, "rb");
        header = sam_hdr_read(bam);
        idx = bam_index_load(bam_fname);

        // DB setup
        // possibly add in multithreading flags
        int sqlite_code = sqlite3_open(dest_fname, &conn);
        if(sqlite_code){
            fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(conn));
            exit(0);
        } else {
            fprintf(stdout, "Opened database successfully\n");
        }
        create_align_table();
        char* err_msg = 0;
        //turn off synchronous writing to disk for increased insertion speed
        sqlite3_exec(conn, "PRAGMA synchronous = OFF", NULL, NULL, &err_msg);

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

// Read in barcodes file and load sequences into barcodes vector
void BamDB::set_barcodes(const char* fname, vector<vector<int> >& vec_p) {
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

int BamDB::create_align_table() {
     char* err_msg = 0;
     int rc;
     //char statement[BUFFER_SIZE];
     char* read1_sql = sqlite3_mprintf("CREATE TABLE IF NOT EXISTS read1 (instrument text, flowcell text, cluster text, tid int, hpos int, tpos int, strand int, bc int, umi int);");
     char* read2_sql = sqlite3_mprintf("CREATE TABLE IF NOT EXISTS read2 (instrument text, flowcell text, cluster text, tid int, hpos int, tpos int, strand int);");


    rc = sqlite3_exec(conn, read1_sql, NULL, NULL, &err_msg);
    if( rc != SQLITE_OK ){
       fprintf(stderr, "SQL error (%d) on read1 table creation: %s\n", rc, err_msg);
    } else{
       fprintf(stdout, "Read1 table created successfully\n");
    }
    sqlite3_free(read1_sql);

    rc = sqlite3_exec(conn, read2_sql, NULL, NULL, &err_msg);
    if( rc != SQLITE_OK ){
       fprintf(stderr, "SQL error (%d) on read2 table creation: %s\n", rc, err_msg);
    } else{
       fprintf(stdout, "Read2 table created successfully\n");
    }
    sqlite3_free(read2_sql);
    return 0;
}

void BamDB::create_reftable() {
    char* err_msg = 0;
    const char* tail = 0;
    int rc;
    char ** names = header->target_name;

    const char* create_table = "CREATE TABLE IF NOT EXISTS reference (name text, tid int);";
    rc = sqlite3_exec(conn, create_table, NULL, NULL, &err_msg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error (%d) on reftable creation: %s \n", rc, err_msg);
    } else {
        cout << "Reference table created successfully." << endl;
    }

    sqlite3_stmt* ref_stmt;
    char sql[BUFFER_SIZE];
    sprintf(sql, "INSERT INTO reference (name, tid) VALUES (?, ?);");
    if (sqlite3_prepare_v2(conn, sql, BUFFER_SIZE, &ref_stmt, &tail) == SQLITE_OK) {
        for (int i = 0; i < header->n_targets; ++i) {
            sqlite3_bind_text(ref_stmt, 1, names[i], -1, SQLITE_TRANSIENT);
            sqlite3_bind_int(ref_stmt, 2, i);
            sqlite3_step(ref_stmt);
            sqlite3_reset(ref_stmt);
        }
        sqlite3_finalize(ref_stmt);
    } else {
        fprintf(stderr, "SQL error on reftable statement prep.\n");
    }

}

void BamDB::index_cluster() {
    char* err_msg = 0;
    int result = 0;

    cout << "Indexing cluster.." << endl;
    if ((result = sqlite3_exec(conn, "CREATE UNIQUE INDEX cluster_index1 ON read1 (cluster)", NULL, NULL, &err_msg)) != SQLITE_OK) {
        fprintf(stderr, "SQL error (%d) on read1.cluster_index: %s \n", result, err_msg);
    }
    if ((result = sqlite3_exec(conn, "CREATE UNIQUE INDEX cluster_index2 ON read2 (cluster)", NULL, NULL, &err_msg)) != SQLITE_OK) {
        fprintf(stderr, "SQL error (%d) on read2.cluster_index: %s \n", result, err_msg);
    }
    sqlite3_free(err_msg);
}

int create_index(BamDB* bamdb) {
    char* err_msg = 0;
    sqlite3_exec(bamdb->get_conn(), "CREATE INDEX name ON reference(name);", NULL, NULL, &err_msg);
    sqlite3_exec(bamdb->get_conn(), "CREATE INDEX bc_hpos ON align(bc, tid, hpos);", NULL, NULL, &err_msg);
    return 0;
}

int process_read1(BamDB* bamdb, dbRecord1* record ,bam1_t* b) {
    //DEBUG_LOG("process read1");
    // continue if read has indels or skipped references
    if (bad_cigar(b)) return 2;

    // continue if mapped to more than one position
    if (filter_multi_reads(b)) return 3;
    record->split_qname(b);
    record->set_positions(b);

    int used_offset;
    record->set_bc(bamdb, b, &used_offset);
    record->set_umi(bamdb, b, used_offset);

    record->insert_to_db();
    return 0;
}

void process_read2(BamDB* bamdb, dbRecord2* record, bam1_t* b) {
    //DEBUG_LOG("process read2");
    record->split_qname(b);
    record->set_positions(b);
    record->insert_to_db();
}

int fill_reads_tid(BamDB* bamdb, dbRecord1* record1, dbRecord2* record2, hts_itr_t* bam_itr) {
    bam1_t* b = bam_init1();
    int result;
    char* err_msg = 0;

    // dbRecord1* record1 = new dbRecord1(bamdb);
    // dbRecord2* record2 = new dbRecord2(bamdb);
    // record1->set_tid(tid)
    sqlite3_exec(bamdb->get_conn(), "BEGIN TRANSACTION", NULL, NULL, &err_msg);

    while ((result = sam_itr_next(bamdb->get_bam(), bam_itr, b)) >= 0) {
        if (b->core.flag&BAM_FMUNMAP) continue;
        if ((b->core.flag&BAM_FREAD1) != 0) {
            process_read1(bamdb, record1, b);
        } else {
            process_read2(bamdb, record2, b);
        }
    }
    bam_destroy1(b);
    sqlite3_exec(bamdb->get_conn(), "END TRANSACTION", NULL, NULL, &err_msg);
    return 0;
}

void fill_reads(BamDB* bamdb) {
    hts_itr_t* bam_itr;
    for (int tid = 0; tid < bamdb->get_header()->n_targets; ++tid) {
        bam_itr = bam_itr_queryi(bamdb->get_idx(), tid, 0, bamdb->get_rlen(tid));

        dbRecord1* record1 = new dbRecord1(bamdb);
        dbRecord2* record2 = new dbRecord2(bamdb);
        record1->set_tid(tid);
        record2->set_tid(tid);

        fill_reads_tid(bamdb, record1, record2, bam_itr);
        delete record1;
        delete record2;
        bam_itr_destroy(bam_itr);
    }
}

void merge_reads(BamDB* bamdb) {
    return;
}

int fill_db(BamDB* bamdb) {
    fill_reads(bamdb);
    bamdb->index_cluster();
    merge_reads(bamdb);
    return 0;
}
