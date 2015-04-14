/***************************
*
*    Table organization
*    merge - Merged reads, full info
         instrument text
         flowcell text
         cluster text
         chrom text
         tid int
         hpos1 int
         tpos1 int
         hpos2 int
         tpos2 int
         isize int
         strand int
         bc int
         umi int
         INDICES:

     collapsed - Reads collapsed to positions defined by
                 read1 and read2 heads
         bc int
         tid int
         lpos1 int
         rpos2 int
         isize int
         strand int
         pos_id int -- tmp.rowid

     pos_rtree - R*Tree for collapsed tid, pos, and strand
         id -- correspondds to collapsed rowid
         tid1
         tid2
         lpos1
         rpos2
         strand1
         strand2

     tmp -
         bc int
         tid int
         hpos1 int
         strand int
         INDICES:
             bc, tid, hpos1, strand

     grouped - UMI counts by reaad1 pos
         //tid int
         //hpos1 int
         //strand int
         pos_id int -- from tmp.rowid
         unique_umi int
         total_umi int\
         INDICES:
             pos_id

General approach:
    Find collapsed positions that intersect with annotation gene set.
    Identify grouped rows that correspond to identified collapsed positions
    Filter grouped rows with at least X number of supporting positions
    Sum unique_umi of grouped rows
*/


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
        //execute(&conns[0], "PRAGMA journal_mode = MEMORY");
        //execute(&conns[0], "PRAGMA synchronous = OFF");
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

void BamDB::create_align_table() {
     char* read1_sql = sqlite3_mprintf("CREATE TABLE IF NOT EXISTS read1 (instrument text, flowcell text, cluster text, chrom text, tid int, hpos int, tpos int, strand int, bc int, umi int);");
     char* read2_sql = sqlite3_mprintf("CREATE TABLE IF NOT EXISTS read2 (instrument text, flowcell text, cluster text, chrom text, tid int, hpos int, tpos int, strand int, isize int);");

    try {
        //sqlite3_exec(conn, read1_sql, NULL, NULL, NULL);
        execute(conns[0], read1_sql);
        execute(conns[0], read2_sql);
    } catch (sql_exception &e) {
        throw e;
        //cerr << "Read table creation error, " <<  e.what() << endl;
        //exit(1);
    }
    cerr << "Read tables created." << endl;
}

void BamDB::create_reftable() {
    int result = 0;
    const char* tail = 0;
    int rc;
    char ** names = header->target_name;

    const char* reftable_sql = "CREATE TABLE IF NOT EXISTS reference (name text, tid int);";
    execute(conns[1], reftable_sql);
    cerr << "Reference table created successfully." << endl;

    sqlite3_stmt* stmt;
    char sql[BUFFER_SIZE];
    sprintf(sql, "INSERT INTO reference (name, tid) VALUES (?, ?);");

        sqlite3_prepare_v2(conns[1], sql, BUFFER_SIZE, &stmt, NULL);
        for (int i = 0; i < header->n_targets; ++i) {
            sqlite3_bind_text(stmt, 1, names[i], -1, SQLITE_TRANSIENT);
            sqlite3_bind_int(stmt, 2, i);
            if ((result = sqlite3_step(stmt)) < 100) {
                throw sql_exception( result, sqlite3_errmsg(conns[0]));
            };
            sqlite3_reset(stmt);
        }
        sqlite3_finalize(stmt);
}
void BamDB::remove_tmp_files() {
    const boost::filesystem::path p(tmp_path);
    if (!boost::filesystem::remove(p)) {
        cerr << "! Failed to clean up temp file: " << p << endl;
    }
}

void BamDB::close_conn(int index) {
    int result = 0;
    if ((result = sqlite3_close(conns[index])) != SQLITE_OK) {
        throw sql_exception(result, sqlite3_errmsg(conns[index]), "error closing");
    }
}

void BamDB::index_cluster() {
    char* err_msg = 0;
    int result = 0;

    cout << "Indexing cluster.." << endl;
    if ((result = sqlite3_exec(conns[0], "CREATE INDEX cluster_index1 ON read1 (cluster)", NULL, NULL, &err_msg)) != SQLITE_OK) {
        fprintf(stderr, "SQL error (%d) on read1.cluster_index: %s \n", result, err_msg);
    }
    if ((result = sqlite3_exec(conns[0], "CREATE INDEX cluster_index2 ON read2 (cluster)", NULL, NULL, &err_msg)) != SQLITE_OK) {
        fprintf(stderr, "SQL error (%d) on read2.cluster_index: %s \n", result, err_msg);
    }
    sqlite3_free(err_msg);
}

void BamDB::merge_tables() {
    cout << "Merging reads..." << endl;

    int result = 0;
    string attach_db = "ATTACH DATABASE '%s' AS merge_db;";
    string attach_db_form = str(boost::format(attach_db) % merge_path.string());
    const char* attach_db_c = attach_db_form.c_str();
    if ((result = sqlite3_exec(conns[0], attach_db_c, NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conns[0]), "attaching db"};
    }

    const char* create_table = "CREATE TABLE IF NOT EXISTS merge_db.merge\
                                    (instrument text,\
                                    flowcell text,\
                                    cluster text,\
                                    chrom text,\
                                    tid int,\
                                    hpos1 int,\
                                    tpos1 int,\
                                    hpos2 int,\
                                    tpos2 int,\
                                    isize int,\
                                    strand int,\
                                    bc int,\
                                    umi int);";
    execute(conns[0], create_table);
    // if ((result = sqlite3_exec(conns[0], create_table, NULL, NULL, NULL)) != SQLITE_OK) {
    //     throw sql_exception { result, sqlite3_errmsg(conns[1]), "merge table create"};
    // }
    sqlite3_exec(conns[0], "BEGIN", NULL, NULL, NULL);
    const char* insert_merge = "INSERT INTO merge_db.merge\
                                SELECT\
                                    read1.instrument,\
                                    read1.flowcell,\
                                    read1.cluster,\
                                    read1.chrom,\
                                    read1.tid,\
                                    read1.hpos,\
                                    read1.tpos,\
                                    read2.hpos,\
                                    read2.tpos,\
                                    read2.isize,\
                                    read1.strand,\
                                    read1.bc,\
                                    read1.umi\
                               FROM read1 JOIN read2 ON read1.cluster=read2.cluster;";
    if ((result = sqlite3_exec(conns[0], insert_merge, NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conns[1]), "insert merge"};
    }
    sqlite3_exec(conns[0], "COMMIT", NULL, NULL, NULL);

    cout << "Merge table created." << endl;
    cout << "Closing tmp database." << endl;
    sqlite3_close(conns[0]);
}

void BamDB::index_merge() {
    open_connection(merge_path.string().c_str(), &conns[1], true);
    //execute(conns[1], "CREATE INDEX merge_positions ON merge(lpos2, rpos1)");
    //execute(conns[1], "CREATE INDEX merge_isize ON merge(isize)");
    execute(conns[1], "CREATE INDEX merge_full ON merge(hpos1, chrom, bc, strand)");
}

void BamDB::collapse_positions() {

    cout << "Collapsing reads to 'collapsed'" << endl;
    int result = 0;

    open_connection(merge_path.string().c_str(), &conns[1], true);
    execute(conns[1], "BEGIN");
    execute(conns[1], "DROP TABLE IF EXISTS collapsed;");
    // sqlite3_exec(conns[1], "BEGIN", NULL, NULL, NULL);
    // const char* drop_table = "DROP TABLE IF EXISTS collapsed;";
    // if ((result = sqlite3_exec(conns[1], drop_table, NULL, NULL, NULL)) != SQLITE_OK) {
    //     throw sql_exception { result, sqlite3_errmsg(conns[1])};

    // }

    const char* create_table = "CREATE TABLE IF NOT EXISTS collapsed (\
                                        bc int,\
                                        tid int,\
                                        lpos1 int,\
                                        rpos2 int,\
                                        isize int,\
                                        strand int,\
                                        pos_id int);";
    execute(conns[1], create_table);
    // if ((result = sqlite3_exec(conns[1], create_table, NULL, NULL, NULL)) != SQLITE_OK) {
    //    // DEBUG_LOG("Create table error");
    //     throw sql_exception { result, sqlite3_errmsg(conns[1])};

    // }
    //sqlite3_open(dest_merge_c, &conns[1])
//const char* sql = read_sql("/sql/aggregate_umi3.sql");


    /* Group by bc, tid, hpos1, strand */
    execute(conns[1], "CREATE TEMP TABLE tmp AS SELECT bc, tid, hpos1, strand FROM merge\
                       GROUP BY\
                           bc,\
                           tid,\
                           hpos1, \
                           strand;");

    execute(conns[1], "CREATE INDEX tmp_pos ON tmp(bc,tid,hpos1,strand);");

    const char* insert_sql = "INSERT INTO collapsed\
                       SELECT\
                           merge.bc,\
                           merge.tid,\
                           CASE WHEN merge.strand IS 0 THEN merge.hpos1 ELSE merge.hpos2 END as lpos1,\
                           CASE WHEN merge.strand IS 0 THEN merge.hpos2 ELSE merge.hpos1 END as rpos2,\
                           isize, \
                           merge.strand,\
                           tmp.rowid\
                       FROM\
                           merge JOIN tmp\
                           WHERE merge.bc = tmp.bc \
                                 AND merge.tid = tmp.tid \
                                 AND merge.hpos1 = tmp.hpos1 \
                                 AND merge.strand = tmp.strand \
                       GROUP BY       \
                           merge.bc,\
                           merge.tid,\
                           merge.hpos1,\
                           merge.hpos2,\
                           merge.strand;";
    execute(conns[1], insert_sql);

    // const char* pos_id = "UPDATE collapsed SET COLUMN pos_id = (SELECT\
    //                           tmp.rowid FROM tmp WHERE\
    //                           bc=bc,\
    //                           tid=tid,\
    //                           CASE


    //"
    // if ((result = sqlite3_exec(conns[1], insert_sql, NULL, NULL, NULL)) != SQLITE_OK) {
    //    // DEBUG_LOG("exec error");
    //     throw sql_exception { result, sqlite3_errmsg(conns[1]), "insertion to collapsed"};

    // }
    sqlite3_exec(conns[1], "COMMIT", NULL, NULL, NULL);
}

void BamDB::create_pos_indices() {
    int result = 0;
    cout << "Indexing positions..." << endl;
    execute(conns[1], "CREATE INDEX pos_index1 ON merge (bc, chrom, hpos1, hpos2, strand)");
    execute(conns[1], "CREATE INDEX pos_index2 ON collapsed (bc, chrom, lpos1, rpos2, strand)");
}

void BamDB::create_idcollapsed() {
    int result = 0;
    cout << "Adding idcollapse to merge..." << endl;
    execute(conns[1], "ALTER TABLE merge ADD COLUMN idcollapsed int");
    // if ((result = sqlite3_exec(conns[1], "ALTER TABLE merge ADD COLUMN idcollapsed int", NULL, NULL, NULL)) != SQLITE_OK) {
    //     throw sql_exception { result, sqlite3_errmsg(conns[1]), "add idcollapsed" };
    // }
    const char* update = "UPDATE merge SET idcollapsed = \
                              (SELECT rowid FROM collapsed WHERE\
                                  merge.bc=collapsed.bc AND\
                                  merge.chrom=collapsed.chrom AND\
                                  merge.hpos1=collapsed.lpos1 AND\
                                  merge.hpos2=collapsed.rpos2 AND\
                                  merge.strand=collapsed.strand);";
    execute(conns[1], update);
    // if ((result = sqlite3_exec(conns[1], update, NULL, NULL, NULL)) != SQLITE_OK) {
    //     throw sql_exception { result, sqlite3_errmsg(conns[1]), "update idcollapsed" };
    // }
}

void BamDB::create_rtree() {
    cerr << "Creating rtree..." << endl;
    execute(conns[1], "CREATE VIRTUAL TABLE pos_rtree USING rtree_i32(id, tid1, tid2, lpos1, rpos2, strand1, strand2);");
    execute(conns[1], "BEGIN TRANSACTION");
    execute(conns[1], "INSERT INTO pos_rtree SELECT rowid, tid, tid, lpos1, rpos2, strand, strand FROM collapsed;");
    execute(conns[1], "END TRANSACTION");
    cerr << "Finished creating rtree..." << endl;
}

void BamDB::create_collapsed_index() {
    cout << "Creating indices..." << endl;
    try {
        //execute(conns[1], "CREATE INDEX chrom_index ON collapsed(chrom);" );
        //execute(conns[1], "CREATE INDEX collapsed_index ON collapsed(bc, chrom, lpos1, lpos2, strand);" );
        execute(conns[1], "CREATE INDEX collapsed_isize ON collapsed(isize)");
        execute(conns[1], "CREATE INDEX collapsed_index ON collapsed(lpos1, rpos2, chrom, lpos1, lpos2, strand);" );
    } catch (sql_exception &e) {
        throw e;
    }
}

void BamDB::group_umi() {
    cout << "Grouping umi..." << endl;

    // const char* create_table = "CREATE TABLE grouped (\
    //                                 tid int,\
    //                                 hpos1 int,\
    //                                 strand int,\
    //                                 unique_umi int,\
    //                                 total_umi int\
    //                             )";

    const char* create_table = "CREATE TABLE grouped (\
                                    pos_id int,\
                                    unique_umi int,\
                                    total_umi int,\
                                    FOREIGN KEY(pos_id)\
                                        REFERENCES collapsed(pos_id)\
                                )";

    execute(conns[1], create_table);

    const char* group_umi = "INSERT INTO grouped\
                             SELECT\
                                tmp.rowid as pos_id ,\
                                count(distinct umi) as unique_umi,\
                                count(umi) as total_umi\
                                FROM\
                                    merge JOIN tmp ON\
                                        merge.bc = tmp.bc\
                                        AND merge.tid = tmp.tid\
                                        AND merge.hpos1 = tmp.hpos1\
                                        AND merge.strand = tmp.strand        \
                                GROUP BY\
                                    tmp.rowid\
                                ";

    execute(conns[1], "BEGIN TRANSACTION");
    execute(conns[1], group_umi);
    execute(conns[1], "END TRANSACTION");
    execute(conns[1], "CREATE INDEX grouped_pos ON grouped(pos_id);");
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
    record->insert_to_db();

    return 0;
}

int process_read2(BamDB* bamdb, dbRecord2* record, bam1_t* b) {
   // DEBUG_LOG("process read2");
   // continue if read has indels or skipped references
    if (bad_cigar(b)) return 2;

    // continue if mapped to more than one position
    if (filter_multi_reads(b)) return 3;
    record->split_qname(b);
    record->set_positions(b);
    record->set_isize(abs(b->core.isize));
    record->insert_to_db();
}

int fill_reads_tid(BamDB* bamdb, dbRecord1* record1, dbRecord2* record2, hts_itr_t* bam_itr) {
   // DEBUG_LOG("fill reads tid");
    bam1_t* b = bam_init1();
    int result;
    //char* err_msg = 0;

    //sqlite3_exec(bamdb->get_conn(), "BEGIN TRANSACTION", NULL, NULL, &err_msg);
    execute(bamdb->get_conn(), "BEGIN");
    while ((result = sam_itr_next(bamdb->get_bam(), bam_itr, b)) >= 0) {
        if (b->core.flag&BAM_FMUNMAP) continue;
        if ((b->core.flag&BAM_FREAD1) != 0) {
            process_read1(bamdb, record1, b);
        } else {
            process_read2(bamdb, record2, b);
        }
    }
    bam_destroy1(b);
    execute(bamdb->get_conn(), "COMMIT");
    //sqlite3_exec(bamdb->get_conn(), "END TRANSACTION", NULL, NULL, &err_msg);
    //sqlite3_free(err_msg);
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

    // bamdb->index_merge();
    // print_time(start);
    //bamdb->close_conn(0);
    bamdb->collapse_positions();
    print_time(start);

    bamdb->create_rtree();
    print_time(start);

    bamdb->group_umi();
    print_time(start);

    bamdb->create_reftable();
    print_time(start);
    // bamdb->create_idcollapsed();
    // print_time(start);

    // bamdb->create_pos_indices();
    // print_time(start);

    // bamdb->create_collapsed_index();
    // print_time(start);

    const boost::filesystem::path m(bamdb->get_merge_path());
    const boost::filesystem::path d(bamdb->get_dest_path());
    boost::filesystem::rename(m, d);

   }   catch (sql_exception &e) {
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
