#include <bam_sql_merge.h>

using namespace std;

// Main BamDB constructor
BamDB::BamDB(const char* bam_fname, const char* dest_fname, const char* barcodes_fname, int umi_length, int bc_min_qual, bool paired_end)
    : bam_fname(bam_fname)
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

    boost::filesystem::path dest_fname_path(dest_fname);
    dest_path = boost::filesystem::complete(dest_fname_path);

    tmp_path = dest_path;
    tmp_path += ".tmp";

    merge_path = dest_path.parent_path();
    merge_path /= "merge.db";

    try {
        open_connection(tmp_path.string().c_str(), &conns[0], true);
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

BamDB::~BamDB() {
    sqlite3_close(conns[0]);
    sqlite3_close(conns[1]);
    sam_close(bam);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    remove_tmp_files();
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
        //exit(1);
    }
    cerr << "Read tables created." << endl;
    return 0;
}

void BamDB::remove_tmp_files() {
    const boost::filesystem::path p(tmp_path);
    if (!boost::filesystem::remove(p)) {
        cerr << "! Failed to clean up temp file: " << p << endl;
    }
}

void BamDB::close_conn(int index) {
    try {
        sqlite3_close(conns[index]);
    } catch (sql_exception &e) {
        cerr << printf("Error closing connection %d, ", index) << e.what() << endl;
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
    cout << "Merging reads..." << endl;

    int result = 0;
    string attach_db = "ATTACH DATABASE '%s' AS merge_db;";
    string attach_db_form = str(boost::format(attach_db) % merge_path.string());
    const char* attach_db_c = attach_db_form.c_str();
    DEBUG_LOG(attach_db_c);
    if ((result = sqlite3_exec(conns[0], attach_db_c, NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conns[0]), "attaching db"};
    }

    const char* create_table = "CREATE TABLE IF NOT EXISTS merge_db.merge (instrument text, flowcell text, cluster text, chrom text, tid int, lpos1 int, lpos2 int, rpos1, rpos2, strand int, bc int, umi int);";
    if ((result = sqlite3_exec(conns[0], create_table, NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conns[1]), "merge table create"};
    }
    sqlite3_exec(conns[0], "BEGIN", NULL, NULL, NULL);
    const char* insert_merge = "INSERT INTO merge_db.merge\
                                SELECT\
                                    read1.instrument,\
                                    read1.flowcell,\
                                    read1.cluster,\
                                    read1.chrom,\
                                    read1.tid,\
                                    CASE WHEN read1.strand IS 0 THEN read1.hpos ELSE read2.hpos END,\
                                    CASE WHEN read1.strand IS 0 THEN read1.tpos ELSE read2.tpos END,\
                                    CASE WHEN read1.strand IS 0 THEN read2.tpos ELSE read1.tpos END,\
                                    CASE WHEN read1.strand IS 0 THEN read2.hpos ELSE read1.hpos END,\
                                    read1.strand,\
                                    read1.bc,\
                                    read1.umi\
                               FROM read1 LEFT JOIN read2 ON read1.cluster=read2.cluster;";
    if ((result = sqlite3_exec(conns[0], insert_merge, NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conns[1]), "insert merge"};
    }

        sqlite3_exec(conns[0], "COMMIT", NULL, NULL, NULL);

    cout << "Merge table created." << endl;
    cout << "Closing tmp database." << endl;
    sqlite3_close(conns[0]);

    return 1;
}

void BamDB::aggregate_umi() {

    cout << "Aggregating reads to 'collapsed'" << endl;
    int result = 0;

    if ((result = sqlite3_open(merge_path.string().c_str(), &conns[1])) != SQLITE_OK) {
        cerr << "! Error opening merge.db" << endl;
        throw sql_exception { result, sqlite3_errmsg(conns[1])};

    }

    sqlite3_exec(conns[1], "BEGIN", NULL, NULL, NULL);
    const char* drop_table = "DROP TABLE IF EXISTS collapsed;";
    if ((result = sqlite3_exec(conns[1], drop_table, NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conns[1])};

    }
    const char* create_table = "CREATE TABLE IF NOT EXISTS collapsed (\
                                        bc int,\
                                        chrom text,\
                                        lpos1 int,\
                                        lpos2 int,\
                                        rpos1 int,\
                                        rpos2 int,\
                                        strand int,\
                                        total_umi int);";
    if ((result = sqlite3_exec(conns[1], create_table, NULL, NULL, NULL)) != SQLITE_OK) {
       // DEBUG_LOG("Create table error");
        throw sql_exception { result, sqlite3_errmsg(conns[1])};

    }
    //sqlite3_open(dest_merge_c, &conns[1])
//const char* sql = read_sql("/sql/aggregate_umi3.sql");
    const char* sql = "INSERT INTO collapsed\
                       SELECT\
                           bc,\
                           chrom,\
                           lpos1,\
                           lpos2,\
                           rpos1,\
                           rpos2,\
                           strand,\
                           count(distinct umi)\
                       FROM\
                           merge\
                       GROUP BY\
                           bc,\
                           chrom,\
                           lpos1,\
                           rpos2,\
                           strand;";

    if ((result = sqlite3_exec(conns[1], sql, NULL, NULL, NULL)) != SQLITE_OK) {
       // DEBUG_LOG("exec error");
        throw sql_exception { result, sqlite3_errmsg(conns[1]), "insertion to collapsed"};

    }
    sqlite3_exec(conns[1], "COMMIT", NULL, NULL, NULL);
}

int BamDB::create_pos_indices() {
    int result = 0;
    cout << "Indexing positions..." << endl;
    if ((result = sqlite3_exec(conns[1], "CREATE INDEX pos_index1 ON merge (bc, chrom, lpos1, rpos2, strand)", NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conns[1]), "position index1" };
    }
    if ((result = sqlite3_exec(conns[1], "CREATE INDEX pos_index2 ON collapsed (bc, chrom, lpos1, rpos2, strand)", NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conns[1]), "position index2" };
    }
    return 0;
}

int BamDB::create_idcollapsed() {
    int result = 0;
    cout << "Adding idcollapse to merge..." << endl;
    if ((result = sqlite3_exec(conns[1], "ALTER TABLE merge ADD COLUMN idcollapsed int", NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conns[1]), "add idcollapsed" };
    }
    const char* update = "UPDATE merge SET idcollapsed = \
                              (SELECT rowid FROM collapsed WHERE\
                                  merge.bc=collapsed.bc AND\
                                  merge.chrom=collapsed.chrom AND\
                                  merge.lpos1=collapsed.lpos1 AND\
                                  merge.rpos2=collapsed.rpos2 AND\
                                  merge.strand=collapsed.strand);";
    if ((result = sqlite3_exec(conns[1], update, NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conns[1]), "update idcollapsed" };
    }
    return 0;
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
   // DEBUG_LOG("process read1");
    // continue if read has indels or skipped references
    if (bad_cigar(b)) return 2;

    // continue if mapped to more than one position
    if (filter_multi_reads(b)) return 3;
    record->split_qname(b);
    record->set_positions(b);

    int used_offset = 0;
    record->set_bc(bamdb, b, &used_offset);
    record->set_umi(bamdb, b, used_offset);

    try {
        record->insert_to_db();
    } catch (sql_exception &e) {
        cerr << "Read1 insertion error, " << e.what() << endl;

    }
    return 0;
}

void process_read2(BamDB* bamdb, dbRecord2* record, bam1_t* b) {
   // DEBUG_LOG("process read2");
    record->split_qname(b);
    record->set_positions(b);

    try {
        record->insert_to_db();
    } catch (sql_exception &e) {
        cerr << "Read2 insertion error, " << e.what() << endl;

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
   //// DEBUG_LOG("fill reads");
   cout << "Filling read tables..." << endl;
    hts_itr_t* bam_itr;
    for (int tid = 0; tid < bamdb->get_header()->n_targets; ++tid) {
        bam_itr = bam_itr_queryi(bamdb->get_idx(), tid, 0, bamdb->get_rlen(tid));
        //DEBUG_LOG(tid);
        dbRecord1* record1 = new dbRecord1(bamdb);
        dbRecord2* record2 = new dbRecord2(bamdb);
       // DEBUG_LOG("Created records");
        record1->set_tid(tid);
        record2->set_tid(tid);
       //// DEBUG_LOG(record1->get_tid());
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

    try {
        fill_reads(bamdb);
        print_time(start);

        bamdb->index_cluster();
        print_time(start);

        bamdb->merge_tables();
        print_time(start);

        bamdb->close_conn(0);
        bamdb->aggregate_umi();
        print_time(start);

        bamdb->create_pos_indices();
        print_time(start);

        bamdb->create_idcollapsed();
        print_time(start);

        bamdb->create_rtree();
        print_time(start);

        bamdb->create_index();
        print_time(start);

        const boost::filesystem::path m(bamdb->get_merge_path());
        const boost::filesystem::path d(bamdb->get_dest_path());
        boost::filesystem::rename(m, d);

    } catch (sql_exception &e) {
        cerr << " ! Error: " << e.what() << endl;
        delete bamdb;
        exit(1);
    } catch (runtime_error &e) {
        cerr << " ! Runtime error: " << e.what() << endl;
        delete bamdb;
        exit(1);
    }
    return 0;
}
