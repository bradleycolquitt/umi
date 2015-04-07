#include <fstream>
#include <sstream>
#include <bam_sql.h>
#include <bam_utils.h>
#include <string_utils.h>
#include <boost/regex.hpp>
#include <gperftools/profiler.h>

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

        // possibly add in multithreading flags
        int sqlite_code = sqlite3_open(dest_fname, &conn);
        if(sqlite_code){
            fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(conn));
            exit(0);
        } else {
            fprintf(stdout, "Opened database successfully\n");
        }

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
        #ifdef DEBUG
        cerr << sequence_pos[2] << " " << sequence_pos[3] << " "
             << bc_length << endl;
        #endif
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

// Create 'reference table' containin human-readable reference names and tid.
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

vector<int> dbRecordSe::insert_to_db() {
    vector<int> codes;
    //codes.resize(12);
    codes.push_back(sqlite3_bind_text(stmt, 1, instrument, -1, SQLITE_TRANSIENT));
    codes.push_back(sqlite3_bind_text(stmt, 2, flowcell, -1, SQLITE_TRANSIENT));
    codes.push_back(sqlite3_bind_text(stmt, 3, cluster, -1, SQLITE_TRANSIENT));
    codes.push_back(sqlite3_bind_int(stmt, 4, tid));
    codes.push_back(sqlite3_bind_int(stmt, 5, read_pos[0]));
    codes.push_back(sqlite3_bind_int(stmt, 6, read_pos[1]));
    codes.push_back(sqlite3_bind_int(stmt, 7, strand));
    codes.push_back(sqlite3_bind_int(stmt, 8, bc));
    codes.push_back(sqlite3_bind_int(stmt, 9, umi));

    codes.push_back(sqlite3_step(stmt));
    //codes.push_back(sqlite3_clear_bindings(stmt));
    codes.push_back(sqlite3_reset(stmt));

    return codes;
}

vector<int> dbRecordPe::insert_to_db() {
    vector<int> kludge;
    return kludge;
}

vector<int> dbRecordPe::insert_to_db(int read_num) {
    vector<int> codes;
    codes.resize(12);
    int result = 0;
    if (read_num == 1) {
        cerr << "insert1" << endl;
        sqlite3_bind_text(insert_stmt, 1, instrument, -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(insert_stmt, 2, flowcell, -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(insert_stmt, 3, cluster, -1, SQLITE_TRANSIENT);
        sqlite3_bind_int(insert_stmt, 4, tid);
        sqlite3_bind_int(insert_stmt, 5, read_pos[0]);
        sqlite3_bind_int(insert_stmt, 6, read_pos[1]);
        sqlite3_bind_int(insert_stmt, 7, 0);
        sqlite3_bind_int(insert_stmt, 8, 0);
        sqlite3_bind_int(insert_stmt, 9, 0);
        sqlite3_bind_int(insert_stmt, 10, strand);
        sqlite3_bind_int(insert_stmt, 11, bc);
        sqlite3_bind_int(insert_stmt, 12, umi);
        #ifdef DEBUG
        cerr << "ins:" << instrument
             << " flow:" << flowcell
             << " cluster:" << cluster
             << " tid:" << tid
             << " head:" << read_pos[0]
             << " tail:" << read_pos[1]
             << " strand:" << strand
             << " bc:" << bc
             << " umi:" << umi << endl;;
        #endif
        if ((result = sqlite3_step(insert_stmt)) != SQLITE_DONE ) {
            fprintf(stderr, "Insertion error (%d): read_num=1, %s\n", result, cluster);
        }
        //sqlite3_clear_bindings(insert_stmt);
        sqlite3_reset(insert_stmt);
    } else if (read_num == 2) {
        cerr << "insert2" << endl;
        sqlite3_bind_text(insert_stmt, 1, "", -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(insert_stmt, 2, "", -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(insert_stmt, 3, "", -1, SQLITE_TRANSIENT);
        sqlite3_bind_int(insert_stmt, 4, 0);
        sqlite3_bind_int(insert_stmt, 5, 0);
        sqlite3_bind_int(insert_stmt, 6, 0);
        sqlite3_bind_int(insert_stmt, 7, read_pos[2]);
        sqlite3_bind_int(insert_stmt, 8, read_pos[3]);
        sqlite3_bind_int(insert_stmt, 9, insert);
        sqlite3_bind_int(insert_stmt, 10, 0);
        sqlite3_bind_int(insert_stmt, 11, 0);
        sqlite3_bind_int(insert_stmt, 12, 0);
        sqlite3_step(insert_stmt);
        //sqlite3_clear_bindings(insert_stmt2);
        sqlite3_reset(insert_stmt);
    }
    return codes;
}

int dbRecordPe::update_record(int read_num) {
    if (read_num == 1) {
        sqlite3_bind_int(update_stmt1, 1, read_pos[0]);
        sqlite3_bind_int(update_stmt1, 2, read_pos[1]);
        sqlite3_bind_int(update_stmt1, 3, strand);
        sqlite3_bind_int(update_stmt1, 4, bc);
        sqlite3_bind_int(update_stmt1, 5, umi);

        sqlite3_step(update_stmt1);
        sqlite3_clear_bindings(update_stmt1);
        sqlite3_reset(update_stmt1);
    } else if (read_num == 2) {
        sqlite3_bind_int(update_stmt2, 1, read_pos[2]);
        sqlite3_bind_int(update_stmt2, 2, read_pos[3]);
        sqlite3_bind_int(update_stmt2, 3, insert);

        sqlite3_step(update_stmt2);
        sqlite3_clear_bindings(update_stmt2);
        sqlite3_reset(update_stmt2);
    }
    return 0;
}

int BamDB::create_align_table() {
     char* err_msg = 0;
     int rc;
     char statement[BUFFER_SIZE];
     if (!paired_end) {
        sprintf(statement, "CREATE TABLE IF NOT EXISTS align (instrument text, flowcell text, cluster text, tid int, hpos int, tpos int, strand int, bc int, umi int);");
     } else {
         sprintf(statement, "CREATE TABLE IF NOT EXISTS align (instrument text, flowcell text, cluster text, tid int, hpos1 int, tpos1 int, hpos2 int, tpos2 int, isize int, strand int, bc int, umi int);");
     }

    rc = sqlite3_exec(conn, statement, NULL, NULL, &err_msg);
    if( rc != SQLITE_OK ){
       fprintf(stderr, "SQL error (%d) on align table creation: %s\n", rc, err_msg);
    } else{
       fprintf(stdout, "Align table created successfully\n");
    }
    sqlite3_free(err_msg);
    return 0;
}

int create_index(BamDB* bamdb) {
    char* err_msg = 0;
    sqlite3_exec(bamdb->get_conn(), "CREATE INDEX name ON reference(name);", NULL, NULL, &err_msg);
    sqlite3_exec(bamdb->get_conn(), "CREATE INDEX bc_hpos ON align(bc, tid, hpos);", NULL, NULL, &err_msg);
    return 0;
}

/***********************
Read Processing utilities
**************************/
int bad_cigar(bam1_t* b) {
    uint32_t* cigar = bam_get_cigar(b);
    int op;
    for (int i=0; i < b->core.n_cigar; ++i) {
        op = bam_cigar_op(cigar[i]);
        if (op == BAM_CINS || op == BAM_CDEL || op == BAM_CREF_SKIP) {
            return 1;
        }
    }
    return 0;
}

// Filter out multimappers (NH>1)
int filter_multi_reads(bam1_t* b) {
    uint8_t *s = bam_aux_get(b, "NH");
    if (s) {
        if (bam_aux2i(s) > 1) {
            //cout << bam_aux2i(s) << endl;
            return 1;
        }
    }
    return 0;
}

/***********************************
Read processing and record updating
************************************/
// Split qname/read metadata
void dbRecord::split_qname(bam1_t* b) {
    char* qname = bam_get_qname(b);
    char* qname_array = strtok(qname, ":");
    int j = 0;
    while (qname_array != NULL) {
        // change these to settors
        if (j == 0) {
            strcpy(instrument, qname_array);
           // record->seinstrument(qname_array);
        } else if (j == 2) {
            strcpy(flowcell, qname_array);
            //record->set_flowcell(qname_array);
        } else if (j== 3) {
            strcpy(cluster, qname_array);
            //record->set_cluster(qname_array);
        } else if (j > 3) {
            strcat(cluster, ":");
            strcat(cluster, qname_array);
        }
        qname_array = strtok(NULL, ":");
        ++j;
    }
}

int compare_barcode_local(vector<vector<int> >::iterator bc_iter, uint8_t* seq, int start, int end) {
    int k = 0;
    int mm = 0;
    for (int j = start; j <= end ; ++j) {
        // #ifdef DEBUG
        // cerr << "j:" << j << " k:" << k << endl;
        // #endif
        if (mm > 1) return 0; // return after second mismatch
        if ((*bc_iter)[k] != bam_seqi(seq, j)) mm += 1;
        ++k;
    }
    return 1; // match found
}

/*
Look for barcode match within given start and end of provided sequence.
Will also search sequence with position offsets provided in bc_offsets vector
*/
int compare_barcode(uint8_t* seq, vector<vector<int> >* barcodes, vector<int>* bc_offsets, int start, int end, int* used_offset) {

    vector<int>::iterator offsets_iter = bc_offsets->begin();
    vector<int>::iterator offsets_iter_end = bc_offsets->end();

    vector<vector<int> >::iterator bc_iter;
    vector<vector<int> >::iterator bc_iter_end = barcodes->end();

    int i = 0;
    for (; offsets_iter != offsets_iter_end; ++offsets_iter) {
        bc_iter = barcodes->begin();
        for (; bc_iter != bc_iter_end ; ++bc_iter) {
            // #ifdef DEBUG
            // cerr << "i:" << i << endl;
            // #endif
            if (compare_barcode_local(bc_iter, seq, start + *offsets_iter, end + *offsets_iter)) {
                *used_offset = *offsets_iter;
                return i;
            }
            ++i;
        }
        i = 0;
    }
    return -1;
}

//for barcodes
int get_sequence(bam1_t* b, int start, int end, vector<vector<int> >* barcodes, \
                 vector<int>* bc_offsets, int min_qual, int* used_offset) {
    uint8_t* seq = bam_get_seq(b);
    uint8_t* qual = bam_get_qual(b);

    // skip read if barcode quality less than min_qual
    int min = 100000;
    int qual_int;
    for (int j = start; j <= end ; ++j) {
        //cout << j << endl;
        qual_int = int(*qual);
        if (qual_int < min) min = qual_int;
        qual += 1;
    }

    if (min < min_qual) return -1;

    // returns index of perfect match or one mismatch
    return compare_barcode(seq, barcodes, bc_offsets, start, end, used_offset);
}

int ShiftAdd(int sum, int digit)
{
     return sum*10 + digit;
}

//intended for umi
int get_sequence(bam1_t* b, int start, int end, int used_offset) {
    uint8_t* seq = bam_get_seq(b);
    uint8_t* qual = bam_get_qual(b);

    start = start + used_offset;
    end = end + used_offset;

    int d = 0;
    // occurrs if sequencing is truncated at beginning
    if (start < 0) {
        start = 0;
        d = abs(used_offset);
    }

    int min = 100000;
    //int qual_int;
    for (int j = start; j <= end ; ++j) {
        if (int(*qual) < min) min = int(*qual);
        qual += 1;
    }

    if (min < 20) {
        return 0;
    } else {
        vector<int> umi_int(end - start + 1 + d, 0);
        // #ifdef DEBUG
        // cerr << "d:" << d << " umi_int size:" << umi_int.size() << endl;
        // #endif
        for (unsigned int j = 0; j < (umi_int.size() - d); ++j) {
            // #ifdef DEBUG
            // cerr << "j:" << j << endl;
            // #endif
             umi_int[j + d] = (int)bam_seqi(seq, j + start);
        }

        // if offset truncates UMI, randomly assign base to truncated positions
        if (d > 0) {
            int nuc_set[] = {1,2,4,8};
            for (int i = 0; i < d; ++i) {
                umi_int[i] = nuc_set[(int)(rand() % (sizeof(nuc_set) / sizeof(int)))];
            }
        }
        return accumulate(umi_int.begin(), umi_int.end(), 0, ShiftAdd);
    }
}

void dbRecord::set_bc(BamDB* bamdb, bam1_t* b, int* used_offset) {
    bc = get_sequence(b, bamdb->get_sequence_pos(2), bamdb->get_sequence_pos(3), bamdb->get_barcodes(), bamdb->get_bc_offsets(), bamdb->get_bc_min_qual(), used_offset);
}

void dbRecord::set_umi(BamDB* bamdb, bam1_t* b, int used_offset) {
    umi = get_sequence(b, bamdb->get_sequence_pos(0), bamdb->get_sequence_pos(1), used_offset);
}

int record_exists(BamDB* bamdb, dbRecordPe* record) {
    //char* err_msg = 0;
    const char* tail = 0;
    char SQL[BUFFER_SIZE];
    sqlite3_stmt* stmt_exists;
    sprintf(SQL, "SELECT EXISTS(SELECT 1 FROM align WHERE cluster=\'%s\'", record->get_cluster());
    sqlite3_prepare_v2(bamdb->get_conn(), SQL, BUFFER_SIZE, &stmt_exists, &tail);
    sqlite3_step(stmt_exists);
    return sqlite3_column_int(stmt_exists, 0);
}

int process_read1(BamDB* bamdb, dbRecord* record ,bam1_t* b) {
    #ifdef DEBUG
    cerr << "process_read1" << endl;
    #endif //DEBUG
    // continue if read has indels or skipped references
    if (bad_cigar(b)) return 2;

    // continue if mapped to more than one position
    if (filter_multi_reads(b)) return 3;

    record->set_positions(b, 1);

    int used_offset;
    record->set_bc(bamdb, b, &used_offset);
    record->set_umi(bamdb, b, used_offset);

    return 0;
}

void process_read2(dbRecordPe* record, bam1_t* b) {
    #ifdef DEBUG
    cerr << "process_read2" << endl;
    #endif //DEBUG
    record->set_positions(b, 2);
    record->set_insert(b);
}

void process_read(BamDB* bamdb, dbRecordSe* record, bam1_t* b) {
    if ((b->core.flag&BAM_FREAD1) != 0) {
        record->split_qname(b);
        process_read1(bamdb, record, b);
        record->insert_to_db();
    }
}


void process_read(BamDB* bamdb, dbRecordPe* record, bam1_t* b) {
    #ifdef DEBUG
    cerr << "process_read" << endl;
    #endif //DEBUG
    record->split_qname(b);
    int rc_exists = record_exists(bamdb, record);
    if ((b->core.flag&BAM_FREAD1) != 0) {
        process_read1(bamdb, record, b);
        if (rc_exists) {
            record->update_record(1);
        } else {
            record->insert_to_db(1);
        }
    } else if ((b->core.flag&BAM_FREAD2) != 0) {
        process_read2(record, b);
        if (rc_exists) {
            record->update_record(2);
            return;
        } else {
            record->insert_to_db(2);
        }
    }

}

int fill_db(BamDB* bamdb) {
    #ifdef DEBUG
    cerr << "fill_db" << endl;
    #endif //DEBUG
    char* err_msg = 0;
    //turn off synchronous writing to disk for increased insertion speed
    sqlite3_exec(bamdb->get_conn(), "PRAGMA synchronous = OFF", NULL, NULL, &err_msg);

    hts_itr_t* bam_itr;
    for (int tid = 0; tid < bamdb->get_header()->n_targets; ++tid) {
        bam_itr = bam_itr_queryi(bamdb->get_idx(), tid, 0, bamdb->get_rlen(tid));
        if (!bamdb->is_pe()) {
            dbRecordSe* record = new dbRecordSe(bamdb);
            record->set_tid(tid);
            fill_db_tid(bamdb, record, bam_itr);
            delete record;
        } else {
            dbRecordPe* record = new dbRecordPe(bamdb);
            record->set_tid(tid);
            fill_db_tid(bamdb, record, bam_itr);
            delete record;
        }
        bam_itr_destroy(bam_itr);
    }
    return 0;
}

template<typename TRecord>
int fill_db_tid(BamDB* bamdb, TRecord* record, hts_itr_t* bam_itr) {
    bam1_t* b = bam_init1();
    int result;
    vector<int> insertion_result;
    char* err_msg = 0;
    sqlite3_exec(bamdb->get_conn(), "BEGIN TRANSACTION", NULL, NULL, &err_msg);

    while ((result = sam_itr_next(bamdb->get_bam(), bam_itr, b)) >= 0) {
        if (b->core.flag&BAM_FMUNMAP) continue;
        process_read(bamdb, record, b);
        // accumulate(insertion_result.begin(), insertion_result.end(),0)
        // if ( > 0) {
        //     cerr << "Insertion errors: ";
        //     for (int i = 0; i < insertion_result.size(); ++i) {
        //         cerr << insertion_result[i] << " ";
        //     }
        //     cerr << endl;
        //     exit(1);
        // }
    }
    bam_destroy1(b);
    sqlite3_exec(bamdb->get_conn(), "END TRANSACTION", NULL, NULL, &err_msg);
    return 0;
}
