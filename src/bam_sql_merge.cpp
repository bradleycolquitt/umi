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

     pos_rtree - R*Tree for collapsed tid, pos, and strand
         id -- correspondds to collapsed rowid
         tid1
         tid2
         lpos1
         rpos2
         strand1
         strand2

     grouped - UMI counts by reaad1 pos
         tid int
         hpos1 int
         strand int
         unique_umi int
         total_umi int
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
BamDB::BamDB(const char* bam_fname, const char* dest_fname, const char* barcodes_fname, int umi_length, int bc_min_qual, bool just_merge, bool just_fill)
    : bam_fname(bam_fname)
    , just_merge(just_merge)
    , just_fill(just_fill)
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

void BamDB::create_align_table() {
     char* read1_sql = sqlite3_mprintf("CREATE TABLE IF NOT EXISTS read1 (instrument text, flowcell text, cluster text, chrom text, tid int, hpos int, tpos int, strand int, bc int, umi int);");
     char* read2_sql = sqlite3_mprintf("CREATE TABLE IF NOT EXISTS read2 (instrument text, flowcell text, cluster text, chrom text, tid int, hpos int, tpos int, strand int, isize int);");

    try {
        execute(conns[0], read1_sql);
        execute(conns[0], read2_sql);
    } catch (sql_exception &e) {
        throw e;
    }
    cerr << "Read tables created." << endl;
}

void BamDB::create_reftable() {
    int result = 0;
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
    execute(conns[1], "CREATE INDEX reference_name ON reference(name);");

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

    close_conn(0);
    open_connection(dest_path.string().c_str(), &conns[1], true);
    const string attach_db = "ATTACH DATABASE '%s' AS read_db;";
    const string attach_db_form = str(boost::format(attach_db) % tmp_path.string());
    const char* attach_db_c = attach_db_form.c_str();
    execute(conns[1], attach_db_c);

    const char* create_table = "CREATE TABLE IF NOT EXISTS merge\
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
    execute(conns[1], create_table);

    execute(conns[1], "BEGIN");
    const char* insert_merge = "INSERT INTO merge\
                                SELECT\
                                    read_db.read1.instrument,\
                                    read_db.read1.flowcell,\
                                    read_db.read1.cluster,\
                                    read_db.read1.chrom,\
                                    read_db.read1.tid,\
                                    read_db.read1.hpos,\
                                    read_db.read1.tpos,\
                                    read_db.read2.hpos,\
                                    read_db.read2.tpos,\
                                    read_db.read2.isize,\
                                    read_db.read1.strand,\
                                    read_db.read1.bc,\
                                    read_db.read1.umi\
                               FROM read_db.read1 JOIN read_db.read2 ON read_db.read1.cluster=read_db.read2.cluster;";

    execute(conns[1], insert_merge);
    execute(conns[1], "COMMIT");
    cout << "Merge table created." << endl;
}

void BamDB::index_merge() {
    execute(conns[1], "CREATE INDEX merge_full ON merge(hpos1, chrom, bc, strand)");
}

void BamDB::collapse_positions() {

    cout << "Collapsing reads to 'collapsed'" << endl;

    execute(conns[1], "BEGIN");
    execute(conns[1], "DROP TABLE IF EXISTS collapsed;");

    const char* create_table = "CREATE TABLE IF NOT EXISTS collapsed (\
                                        bc int,\
                                        tid int,\
                                        lpos1 int,\
                                        rpos2 int,\
                                        isize int,\
                                        strand int,\
                                        pos_id int);";
    execute(conns[1], create_table);

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
                           merge.strand\
                       HAVING isize BETWEEN 150 AND 1000;";
    execute(conns[1], insert_sql);
    execute(conns[1], "COMMIT");
}

void BamDB::create_pos_indices() {
    cout << "Indexing positions..." << endl;
    execute(conns[1], "CREATE INDEX pos_index1 ON merge (bc, chrom, hpos1, hpos2, strand)");
    execute(conns[1], "CREATE INDEX pos_index2 ON collapsed (bc, chrom, lpos1, rpos2, strand)");
}

void BamDB::create_idcollapsed() {
    cout << "Adding idcollapse to merge..." << endl;
    execute(conns[1], "ALTER TABLE merge ADD COLUMN idcollapsed int");

    const char* update = "UPDATE merge SET idcollapsed = \
                              (SELECT rowid FROM collapsed WHERE\
                                  merge.bc=collapsed.bc AND\
                                  merge.chrom=collapsed.chrom AND\
                                  merge.hpos1=collapsed.lpos1 AND\
                                  merge.hpos2=collapsed.rpos2 AND\
                                  merge.strand=collapsed.strand);";
    execute(conns[1], update);
}

void BamDB::create_rtree() {
    cerr << "Creating rtree..." << endl;
    execute(conns[1], "CREATE VIRTUAL TABLE pos_rtree USING rtree_i32(id, tid1, tid2, lpos1, rpos2, strand1, strand2);");
    execute(conns[1], "BEGIN TRANSACTION");
    execute(conns[1], "INSERT INTO pos_rtree SELECT rowid, tid, tid, lpos1, rpos2, strand, strand FROM collapsed WHERE lpos1 < rpos2;");
    execute(conns[1], "END TRANSACTION");
    cerr << "Finished creating rtree..." << endl;
}

void BamDB::create_collapsed_index() {
    cout << "Creating indices..." << endl;
    try {
        execute(conns[1], "CREATE INDEX collapsed_isize ON collapsed(isize)");
        execute(conns[1], "CREATE INDEX collapsed_index ON collapsed(lpos1, rpos2, chrom, lpos1, lpos2, strand);" );
    } catch (sql_exception &e) {
        throw e;
    }
}

void BamDB::group_umi() {
    cout << "Grouping umi..." << endl;

    execute(conns[1], "CREATE TABLE grouped AS SELECT bc, tid, hpos1, strand, count(distinct umi) AS unique_umi, count(umi) AS total_umi FROM merge\
                       GROUP BY\
                           bc,\
                           tid,\
                           hpos1, \
                           strand;");
    execute(conns[1], "CREATE INDEX grouped_pos ON grouped(bc,tid,hpos1,strand);");

    // const char* create_table = "CREATE TABLE grouped (\
    //                                 pos_id int,\
    //                                 unique_umi int,\
    //                                 total_umi int,\
    //                                 FOREIGN KEY(pos_id)\
    //                                     REFERENCES collapsed(pos_id)\
    //                             )";
    // execute(conns[1], create_table);

    // const char* group_umi = "INSERT INTO grouped\
    //                          SELECT\
    //                             tmp.rowid as pos_id ,\
    //                             count(distinct umi) as unique_umi,\
    //                             count(umi) as total_umi\
    //                             FROM\
    //                                 merge JOIN tmp ON\
    //                                     merge.bc = tmp.bc\
    //                                     AND merge.tid = tmp.tid\
    //                                     AND merge.hpos1 = tmp.hpos1\
    //                                     AND merge.strand = tmp.strand        \
    //                             GROUP BY\
    //                                 tmp.rowid\
    //                             ";

    // execute(conns[1], "BEGIN TRANSACTION");
    // execute(conns[1], group_umi);
    // execute(conns[1], "END TRANSACTION");
    // execute(conns[1], "CREATE INDEX grouped_pos ON grouped(pos_id);");
}

int process_read1(BamDB* bamdb, dbRecord1* record ,bam1_t* b) {
    /* continue if read has indels or skipped references */
    if (bad_cigar(b)) return 2;

    /* continue if mapped to more than one position */
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
    /* continue if read has indels or skipped references */
    if (bad_cigar(b)) return 2;

    /* continue if mapped to more than one position */
    if (filter_multi_reads(b)) return 3;
    record->split_qname(b);
    record->set_positions(b);
    record->set_isize(abs(b->core.isize));
    record->insert_to_db();

    return 0;
}

int fill_reads_tid(BamDB* bamdb, dbRecord1* record1, dbRecord2* record2, hts_itr_t* bam_itr) {
    bam1_t* b = bam_init1();
    int result;
    execute(bamdb->get_conn(), "BEGIN");
    while ((result = sam_itr_next(bamdb->get_bam(), bam_itr, b)) >= 0) {
        if (b->core.flag&BAM_FUNMAP) continue;
        if ((b->core.flag&BAM_FREAD1) != 0) {
            process_read1(bamdb, record1, b);
        } else {
            process_read2(bamdb, record2, b);
        }
    }
    bam_destroy1(b);
    execute(bamdb->get_conn(), "COMMIT");
    return 0;
}

void fill_reads(BamDB* bamdb) {
    cout << "Filling read tables..." << endl;
    hts_itr_t* bam_itr;
    for (int tid = 0; tid < bamdb->get_header()->n_targets; ++tid) {
        bam_itr = bam_itr_queryi(bamdb->get_idx(), tid, 0, bamdb->get_rlen(tid));
        dbRecord1* record1 = new dbRecord1(bamdb);
        dbRecord2* record2 = new dbRecord2(bamdb);
        record1->set_tid(tid);
        record2->set_tid(tid);
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
        if (bamdb->fill_only())
        {
            boost::filesystem::
            rename(bamdb->get_tmp_path(), bamdb->get_dest_path());
            return 1;
        }
        bamdb->index_cluster();
        print_time(start);

        bamdb->merge_tables();
        print_time(start);
        if (bamdb->merge_only()) return 1;


     // bamdb->collapse_positions();
        // print_time(start);

        // bamdb->create_rtree();
        // print_time(start);

        bamdb->group_umi();
        print_time(start);

        bamdb->create_reftable();
        print_time(start);

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
