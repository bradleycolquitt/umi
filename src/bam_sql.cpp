#include <fstream>
#include <sstream>
#include <bam_sql.h>
#include <bam_utils.h>
#include <string_utils.h>
#include <boost/regex.hpp>
#include <gperftools/profiler.h>
//#include <seqan/sequence.h>
//#define DEBUG

using namespace std;

BamDB::BamDB(const char* bam_fname, const char* dest_fname, const char* barcodes_fname, int umi_length, int bc_min_qual)
    : bam_fname(bam_fname)
    , dest_fname(dest_fname)
    , bc_min_qual(bc_min_qual)
    {
        bam = sam_open(bam_fname, "rb");
        header = sam_hdr_read(bam);
        idx = bam_index_load(bam_fname);
        total_mapped = count_bam_records(idx, header);

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

        sequence_pos.push_back(0);
        sequence_pos.push_back(umi_length-1);
        sequence_pos.push_back(umi_length);

        size_t bc_length = barcodes[0].size();
        sequence_pos.push_back(sequence_pos[2] + bc_length);

        int offsets_array[] = {0,-1,1}; // defines offsets used during barcode search
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

void BamDB::create_reftable() {
    char* err_msg = 0;
    int rc;
    char ** names = header->target_name;
    rc = sqlite3_exec(conn, "CREATE TABLE IF NOT EXISTS reference (name text, tid int)", NULL, NULL, &err_msg);

    char stmt[BUFFER_SIZE];
    for (int i = 0; i < header->n_targets; ++i) {
        sprintf(stmt, "INSERT INTO reference VALUES (\'%s\', %i)", names[i], i);
        rc = sqlite3_exec(conn, stmt, NULL, NULL, &err_msg);
        if ( rc != SQLITE_OK) {
            fprintf(stderr, "SQL error: %s\n", err_msg);
        }
    }
}

int create_table(BamDB* bamdb) {
    char* err_msg = 0;
    int rc;
    const char* statement = "CREATE TABLE IF NOT EXISTS align (                      \
                                                           instrument text,          \
                                                           flowcell text,            \
                                                           cluster text,             \
                                                           tid int,                  \
                                                           hpos int,                 \
                                                           tpos int,                 \
                                                           strand int,               \
                                                           bc int,                   \
                                                           umi int);";

    rc = sqlite3_exec(bamdb->get_conn(), statement, NULL, NULL, &err_msg);
    if( rc != SQLITE_OK ){
       fprintf(stderr, "SQL error: %s\n", err_msg);
       sqlite3_free(err_msg);
    } else{
       fprintf(stdout, "Align table created successfully\n");
    }
    return 0;
}

int create_index(BamDB* bamdb) {
    char* err_msg = 0;
    sqlite3_exec(bamdb->get_conn(), "CREATE INDEX name ON reference(name);", NULL, NULL, &err_msg);
    sqlite3_exec(bamdb->get_conn(), "CREATE INDEX bc_hpos ON align(bc, tid, hpos);", NULL, NULL, &err_msg);
    return 0;
}

int insert_to_db(dbRecord* record, sqlite3_stmt* stmt) {

    sqlite3_bind_text(stmt, 1, record->instrument, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, record->flowcell, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 3, record->cluster, -1, SQLITE_TRANSIENT);
    sqlite3_bind_int(stmt, 4, record->tid);
    sqlite3_bind_int(stmt, 5, record->pos_head);
    sqlite3_bind_int(stmt, 6, record->pos_tail);
    sqlite3_bind_int(stmt, 7, record->strand);
    sqlite3_bind_int(stmt, 8, record->bc);
    sqlite3_bind_int(stmt, 9, record->umi);

    sqlite3_step(stmt);

    sqlite3_clear_bindings(stmt);
    sqlite3_reset(stmt);

    return 0;
}

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

// Split qname/read metadata
void split_qname(bam1_t* b, dbRecord* record) {
    char* qname = bam_get_qname(b);
    char* qname_array = strtok(qname, ":");
    int j = 0;
    while (qname_array != NULL) {
        if (j == 0) {
            strcpy(record->instrument, qname_array);
        } else if (j == 2) {
            strcpy(record->flowcell, qname_array);
        } else if (j== 3) {
            strcpy(record->cluster, qname_array);
        } else if (j > 3) {
            strcat(record->cluster, ":");
            strcat(record->cluster, qname_array);
        }
        qname_array = strtok(NULL, ":");
        ++j;
    }

}

int compare_barcode_local(vector<vector<int> >::iterator bc_iter, uint8_t* seq, int start, int end) {
    int k = 0;
    int mm = 0;
    for (int j = start; j <= end ; ++j) {
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
int compare_barcode(uint8_t* seq, vector<vector<int> >* barcodes, vector<int>* bc_offsets, int start, int end, int& used_offset) {

    vector<int>::iterator offsets_iter = bc_offsets->begin();
    vector<int>::iterator offsets_iter_end = bc_offsets->end();

    vector<vector<int> >::iterator bc_iter;
    vector<vector<int> >::iterator bc_iter_end = barcodes->end();

    int i = 0;
    for (; offsets_iter != offsets_iter_end; ++offsets_iter) {
        bc_iter = barcodes->begin();
        for (; bc_iter != bc_iter_end ; ++bc_iter) {
            if (compare_barcode_local(bc_iter, seq, start + *offsets_iter, end + *offsets_iter)) {
                used_offset = *offsets_iter;
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
                 vector<int>* bc_offsets, int min_qual, int& used_offset) {
    uint8_t* seq = bam_get_seq(b);
    uint8_t* qual = bam_get_qual(b);

    // skip read if barcode quality less than min_qual
    int min = 100000;
    int qual_int;
    for (int j = start; j <= end ; ++j) {
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
        for (unsigned int j = 0; j < umi_int.size(); ++j) {
             umi_int[j + d] = (int)bam_seqi(seq, j + start);
        }

        // if offset truncates UMI, randomly assign base to truncated positions
        if (d > 0) {
            int nuc_set[] = {1,2,4,8};
            for (int i = 0; i < d; ++i) {
                umi_int[i] = nuc_set[(int)(rand() % (sizeof(nuc_set) / sizeof(int)))];
            }
        }

        // return as vector<int> as int
        return accumulate(umi_int.begin(), umi_int.end(), 0, ShiftAdd);
    }
}

int fill_db(BamDB* bamdb) {
    char* err_msg = 0;

    //turn off synchronous writing to disk for increased insertion speed
    sqlite3_exec(bamdb->get_conn(), "PRAGMA synchronous = OFF", NULL, NULL, &err_msg);

    hts_itr_t* bam_itr;
    for (int tid = 0; tid < bamdb->get_header()->n_targets; ++tid) {
        bam_itr = bam_itr_queryi(bamdb->get_idx(), tid, 0, bamdb->get_rlen(tid));

        fill_db_tid(bamdb, tid, bam_itr);
    }
    sqlite3_close(bamdb->get_conn());
    return 0;

}

int fill_db_tid(BamDB* bamdb, int tid, hts_itr_t* bam_itr) {
    dbRecord record = {};
    record.tid = tid;

    bam1_t* b = bam_init1();
    int result;
    int i = 0;

    sqlite3_stmt* stmt;
    char* err_msg = 0;
    const char* tail = 0;
    char SQL[BUFFER_SIZE];

    vector<int> sequence_pos = bamdb->get_sequence_pos();

    sprintf(SQL, "INSERT INTO align VALUES (@IN, @FL, @CL, @TID, @HPOS, @TPOS, @STR, @BC, @UMI);");
    sqlite3_prepare_v2(bamdb->get_conn(),  SQL, BUFFER_SIZE, &stmt, &tail);

    sqlite3_exec(bamdb->get_conn(), "BEGIN TRANSACTION", NULL, NULL, &err_msg);

    while ((result = sam_itr_next(bamdb->get_bam(), bam_itr, b)) >= 0) {
        if (b->core.flag&BAM_FMUNMAP) continue;
        //bamdb->increment_read();
        // FIX ME:
        // if (bamdb->get_read() % progress_size == 0) {
        //      cout << (float)bamdb->get_read() / bamdb->get_mapped() * 100  << "%..." << flush;
        // }
        //if (i == 100) exit(1);
        //cout << "test" << endl;
        ++i;

        // continue if read has indels or skipped refernences
        if (bad_cigar(b)) continue;

        // continue if mapped to more than one position
        if (filter_multi_reads(b)) {
            continue;
        }

        split_qname(b, &record);

        if (bam_is_rev(b)) {
            record.pos_tail = b->core.pos + 1;
            record.pos_head = bam_endpos(b);
            record.strand = true;
        } else {
            record.pos_head = b->core.pos + 1;
            record.pos_tail = bam_endpos(b);
            record.strand = false;
        }
        // +1 offset for comparison with 1-based indexing

        int used_offset;
        // extract barcode
        record.bc = get_sequence(b, sequence_pos[2], sequence_pos[3], \
                                    bamdb->get_barcodes(), \
                                    bamdb->get_bc_offsets(), \
                                    bamdb->get_bc_min_qual(), \
                                    used_offset);

        // extract umi
        record.umi = get_sequence(b, sequence_pos[0], sequence_pos[1], used_offset);

        try {
            insert_to_db(&record, stmt);
        } catch ( exception& e ) {
            exit(1);
        }
    }

    sqlite3_exec(bamdb->get_conn(), "END TRANSACTION", NULL, NULL, &err_msg);

    return 0;
}
