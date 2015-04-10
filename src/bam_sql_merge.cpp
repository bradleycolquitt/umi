#include <bam_sql_merge.h>

using namespace std;

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

        char ** chrom_names = header->target_name;
        for (int i = 0 ; i < header->n_targets ; ++i) {
            chroms[i] = chrom_names[i];
        }

        // DB setup
        char tmp[BUFFER_SIZE];
        strcpy(tmp, dest_fname);
        strcat(tmp, ".tmp");
        dest_tmp = tmp;

        boost::filesystem::path p(dest_fname);
        boost::filesystem::path dir = p.parent_path();
        string path = dir.string() + "/merge.db";
        char * path_c = path.c_str();
        dest_merge = path_c;

        try {
            open_connection(dest_tmp, &conns[0], false);
        } catch (sql_exception &e) {
            cerr << "! Error: opening databases, " << e.what() << endl;
        }

        create_align_table();

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
     char* read1_sql = sqlite3_mprintf("CREATE TABLE IF NOT EXISTS read1 (instrument text, flowcell text, cluster text, chrom text, tid int, hpos int, tpos int, strand int, bc int, umi int);");
     char* read2_sql = sqlite3_mprintf("CREATE TABLE IF NOT EXISTS read2 (instrument text, flowcell text, cluster text, chrom text, tid int, hpos int, tpos int, strand int);");

    try {
        //sqlite3_exec(conn, read1_sql, NULL, NULL, NULL);
        execute(conns[0], read1_sql);
        execute(conns[0], read2_sql);
    } catch (sql_exception &e) {
        cerr << "Read table creation error, " <<  e.what() << endl;
        exit(1);
    }
    cerr << "Read tables created." << endl;
    return 0;
}

void BamDB::close_conn(int index) {
    try {
        sqlite3_close(conns[index]);
    } catch (sql_exception &e) {
        cerr << printf("Error closing connection %d, ", index) << e.what() << endl;
    }
}
void BamDB::create_reftable() {
    char ** names = header->target_name;

    const char* reftable_sql = "CREATE TABLE IF NOT EXISTS reference (name text, tid int);";
    try {
        execute(conns[1], reftable_sql);
    } catch (sql_exception &e) {
        cerr << "Reference table creation, " << e.what() << endl;
        return;
    }
    cerr << "Reference table created successfully." << endl;

    sqlite3_stmt* stmt;
    char sql[BUFFER_SIZE];
    sprintf(sql, "INSERT INTO reference (name, tid) VALUES (?, ?);");
    try {
        sqlite3_prepare_v2(conns[1], sql, -1, &stmt, NULL);

        //prepare_statment(conn, sql, stmt);
        for (int i = 0; i < header->n_targets; ++i) {
            DEBUG_LOG(names[i]);
            sqlite3_bind_text(stmt, 1, names[i], -1, SQLITE_TRANSIENT);
            //bind(conn, stmt, 1, names[i]);
            bind(conns[1], stmt, 2, i);
            //step(conns[1], stmt);
            sqlite3_clear_bindings(stmt);
            sqlite3_reset(stmt);
        }
        sqlite3_finalize(stmt);
    } catch (sql_exception &e) {
        cerr << "Reference table insertion, " << e.what() << endl;
        return;
    }
}

void BamDB::index_cluster() {
    char* err_msg = 0;
    int result = 0;

    cout << "Indexing cluster.." << endl;
    if ((result = sqlite3_exec(conns[0], "CREATE UNIQUE INDEX cluster_index1 ON read1 (cluster)", NULL, NULL, &err_msg)) != SQLITE_OK) {
        fprintf(stderr, "SQL error (%d) on read1.cluster_index: %s \n", result, err_msg);
    }
    if ((result = sqlite3_exec(conns[0], "CREATE INDEX cluster_index2 ON read2 (cluster)", NULL, NULL, &err_msg)) != SQLITE_OK) {
        fprintf(stderr, "SQL error (%d) on read2.cluster_index: %s \n", result, err_msg);
    }
    sqlite3_free(err_msg);
}

int BamDB::merge_tables() {
    boost::filesystem::path p(get_selfpath());
    boost::filesystem::path dir = p.parent_path();
    string path = dir.string();
    path = path + "/sql/merge_tables2.sql";

    ifstream sql_file(path);
    string sql_contents((istreambuf_iterator<char>(sql_file)), istreambuf_iterator<char>());

    const char* sql = sql_contents.c_str();
    // FIX THIS! const char* sql = read_sql("/sql/merge_tables2.sql");
    //DEBUG_LOG(sql);
    sqlite3_stmt* stmt;
    int result = 0;
    const char* tail;
    try {
        sqlite3_prepare_v2(conns[0], sql, -1, &stmt, &tail);
        sqlite3_bind_text(stmt, 1, dest_fname, -1, SQLITE_TRANSIENT);
//         DEBUG_LOG(stmt);
        result = sqlite3_step(stmt);
        sqlite3_reset(stmt);
        //DEBUG_LOG(tail);
        while (*tail != '\0') {
            //DEBUG_LOG(tail);
            result = sqlite3_prepare_v2(conns[0], tail, -1, &stmt, &tail);
            if (result != SQLITE_OK) {
                throw sql_exception { result, sqlite3_errmsg(conns[0])};
            }
            //DEBUG_LOG(result);
            sqlite3_step(stmt);
            sqlite3_reset(stmt);
        }
        sqlite3_finalize(stmt);

    } catch (sql_exception &e) {
        cerr << "Merge table insertion, " << e.what() << endl;
        sqlite3_finalize(stmt);
        exit(1);

    }

    cerr << "Merge table created." << endl;
    return 1;
}

// void BamDB::drop_read_tables() {
//     cerr << "Dropping read tables..." << endl;
//     try {
//         execute(conn, "DROP TABLE read1;");
//         execute(conn, "DROP TABLE read2;");
//     } catch (sql_exception &e) {
//         cerr << "Error removing read tables, " << e.what() << endl;
//     }
//     cerr << "Finished drop read tables..." << endl;
// }

void BamDB::aggregate_umi() {

    cout << "Aggregating reads to 'collapsed'" << endl;

    open_connection(dest_merge, &conns[1], false);

    const char* sql = read_sql("/sql/aggregate_umi2.sql");
    //cerr << sql << endl;
    sqlite3_stmt* stmt;
    int result = 0;
    const char* tail;

    //try {
    //prepare_statement(conns[1], sql, &stmt, &tail);
    if ((result = sqlite3_prepare_v2(conns[1], sql, -1, &stmt, &tail)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conns[1])};
    }

    if ((result = sqlite3_step(stmt)) <100 ) {
        throw sql_exception { result, sqlite3_errmsg(conns[1])};
    }

    sqlite3_reset(stmt);
    while (*tail != '\0') {
        //DEBUG_LOG(tail);
        try {
        result = sqlite3_prepare_v2(conns[1], tail, -1, &stmt, &tail);
        if (result != SQLITE_OK) {
            throw sql_exception { result, sqlite3_errmsg(conns[1])};
        }
        if ((result = sqlite3_step(stmt)) < 100 ) {
            throw sql_exception { result, sqlite3_errmsg(conns[1])};
        }

        //prepare_statment(&conns[1], tail, &stmt, &tail);
         // if ((result = sqlite3_step(stmt)) < 100 ) {
         //     throw sql_exception { result, sqlite3_errmsg(conns[1])};
         // }
        //step((conns[1]), &stmt);
        reset(&conns[1], &stmt);
        } catch (sql_exception &e) {
            cerr << " ! Error: aggregating, " << e.what() << endl;
        }
        //sqlite3_reset(stmt);
    }
    sqlite3_finalize(stmt);

    // try {
    //     step_multiple(conns[1], sql);
    // } catch (sql_exception &e) {
    //     cerr << " ! Error : Aggregate umi, " << e.what() << endl;
    // }
}

int BamDB::create_rtree() {
    cerr << "Creating rtree..." << endl;
    try {
        execute(conns[1], "CREATE VIRTUAL TABLE pos_rtree USING rtree_i32(id, lpos1, rpos2);");
    } catch (sql_exception &e) {
        cerr << "Error creating, " << e.what() << endl;
    }
    try {
        execute(conns[1], "BEGIN TRANSACTION");
        execute(conns[1], "INSERT INTO pos_rtree SELECT rowid, lpos1, rpos2 FROM collapsed;");
        execute(conns[1], "END TRANSACTION");
    } catch (sql_exception &e) {
        cerr << "Error populating rtree, " << e.what() << endl;
    }
    cerr << "Finished creating rtree..." << endl;

    return 0;
}

int BamDB::create_index() {
    cout << "Creating indices..." << endl;
    try {
        //execute(conns[1], "CREATE INDEX chrom_index ON merge(chrom);" );
        execute(conns[1], "CREATE INDEX chrom_index ON collapsed(chrom);" );
//        execute(conn, "CREATE INDEX info_index ON merge(bc, umi);" );
    } catch (sql_exception &e) {
        cerr << "Error creating indices, " << e.what() << endl;
    }
    //sqlite3_exec(conn, "CREATE INDEX chrom_index ON merge(chrom);", NULL, NULL, NULL);
    //sqlite3_exec(conn, "CREATE INDEX bc_hpos ON align(bc, tid, hpos);", NULL, NULL, &err_msg);
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

    try {
        record->insert_to_db();
    } catch (sql_exception &e) {
        cerr << "Read1 insertion error, " << e.what() << endl;
        exit(1);
    }
    return 0;
}

void process_read2(BamDB* bamdb, dbRecord2* record, bam1_t* b) {
    //DEBUG_LOG("process read2");
    record->split_qname(b);
    record->set_positions(b);

    try {
        record->insert_to_db();
    } catch (sql_exception &e) {
        cerr << "Read2 insertion error, " << e.what() << endl;
        exit(1);
    }
}

int fill_reads_tid(BamDB* bamdb, dbRecord1* record1, dbRecord2* record2, hts_itr_t* bam_itr) {
   // DEBUG_LOG("fill reads tid");
    bam1_t* b = bam_init1();
    int result;
    char* err_msg = 0;

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
    sqlite3_free(err_msg);
    return 0;
}

void fill_reads(BamDB* bamdb) {
   // DEBUG_LOG("fill reads");
   cout << "Filling read tables..." << endl;
    hts_itr_t* bam_itr;
    for (int tid = 0; tid < bamdb->get_header()->n_targets; ++tid) {
        bam_itr = bam_itr_queryi(bamdb->get_idx(), tid, 0, bamdb->get_rlen(tid));
        //DEBUG_LOG(tid);
        dbRecord1* record1 = new dbRecord1(bamdb);
        dbRecord2* record2 = new dbRecord2(bamdb);
        //DEBUG_LOG("Created records");
        record1->set_tid(tid);
        record2->set_tid(tid);
        //DEBUG_LOG(record1->get_tid());
        record1->set_chrom(bamdb, tid);

        record2->set_chrom(bamdb, tid);

        fill_reads_tid(bamdb, record1, record2, bam_itr);
        delete record1;
        delete record2;
        bam_itr_destroy(bam_itr);
    }
}

int fill_db(BamDB* bamdb) {
    clock_t start;
    start = std::clock();

    fill_reads(bamdb);
    print_time(start);
    bamdb->index_cluster();
    print_time(start);
    if (bamdb->merge_tables()) {
        print_time(start);
        bamdb->close_conn(0);
        remove(bamdb->get_tmp_name());
        bamdb->aggregate_umi();
        print_time(start);
        bamdb->create_rtree();
        print_time(start);
        bamdb->create_index();
        print_time(start);
        rename("merge.db", bamdb->get_dest_name());
    } else  {
        cerr << "Merge reads failed!" << endl;
    }
    return 0;
}
