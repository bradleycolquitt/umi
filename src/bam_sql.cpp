#include <fstream>
#include <sstream>
#include <bam_sql.h>
#include <bam_utils.h>
#include <string_utils.h>
#include <boost/regex.hpp>
#include <gperftools/profiler.h>
//#define DEBUG
//#define BUFFER_SIZE 256

using namespace std;

BamDB::BamDB(const char* bam_fname, const char* dest_fname, const char* barcodes_fname)
    : bam_fname(bam_fname)
    , dest_fname(dest_fname)
    {
        bam = sam_open(bam_fname, "rb");
        header = sam_hdr_read(bam);
        idx = bam_index_load(bam_fname);
        total_mapped = count_bam_records(idx, header);

        // // possibly add in multithreading flagB
        int sqlite_code = sqlite3_open(dest_fname, &conn);
        if(sqlite_code){
            fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(conn));
            exit(0);
        } else {
            fprintf(stdout, "Opened database successfully\n");
        }

        // // Barcode setup
        char bc_path[255];
        strcpy(bc_path, BC_PATH);
        strcat(bc_path, barcodes_fname);
        strcat(bc_path, ".txt");
        //cout << bc_path << endl;
        try {
            set_barcodes(bc_path, barcodes);
        } catch (exception &e) {
            cout << e.what() << endl;
        }
    }

// Read in barcodes file and load sequences intobarcodes vector
// consider encoding barcodes as bits or just int
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

        vec_p.push_back(seq2int(sline[1]));
    }
    #ifdef DEBUG
    for (vector<vector<int> >::iterator it = vec_p.begin(); it != vec_p.end(); ++it) {
        for (vector<int>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
        cout << (*it2);
        }
        cout << endl;
    }
    #endif
}

static int callback(void *NotUsed, int argc, char **argv, char **azColName){
   int i;
   for(i=0; i<argc; i++){
       printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
   }
   printf("\n");
   return 0;
}

int create_table(BamDB* bamdb) {
    char* err_msg = 0;
    int rc;
    const char* statement = "CREATE TABLE IF NOT EXISTS align (                      \
                                                           instrument text,          \
                                                           flowcell text,            \
                                                           cluster text,             \
                                                           tid int,                  \
                                                           position int,             \
                                                           strand int,               \
                                                           bc int,                   \
                                                           umi int);";

    rc = sqlite3_exec(bamdb->get_conn(), statement, callback, 0, &err_msg);
    if( rc != SQLITE_OK ){
       fprintf(stderr, "SQL error: %s\n", err_msg);
       sqlite3_free(err_msg);
    } else{
       fprintf(stdout, "Align table created successfully\n");
    }
    return 0;
}

int insert_to_db(BamDB* bamdb, dbRecord* record, sqlite3_stmt* stmt) {

    sqlite3_bind_text(stmt, 1, record->instrument, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, record->flowcell, -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 3, record->cluster, -1, SQLITE_TRANSIENT);
    sqlite3_bind_int(stmt, 4, record->tid);
    sqlite3_bind_int(stmt, 5, record->pos);
    sqlite3_bind_int(stmt, 6, record->strand);
    sqlite3_bind_int(stmt, 7, record->bc);
    sqlite3_bind_int(stmt, 8, record->umi);

    sqlite3_step(stmt);

    sqlite3_clear_bindings(stmt);
    sqlite3_reset(stmt);

    return 0;
}



int bad_cigar(bam1_t* b) {
    uint32_t* cigar = bam_get_cigar(b);
    for (int i=0; i < b->core.n_cigar; ++i) {
        uint32_t op = bam_cigar_op(cigar[i]);
        if (op > 0 & op < 4) {
            return 1;
        }
    }
    return 0;
}

// get aux. filter out NH>1
int filter_bad_reads(bam1_t* b, const char tag[2]) {
        // // FIGURE OUT HOW TO TRANSLATE s
        // uint8_t *s = bam_aux_get(b, "NH");
        // if (s) {
        //     cout << "hit";
        //     char* decoded_nh = bam_aux2Z(s);
        //     cout << decoded_nh;
        //     // for (int i =0; i<10; ++i) {
        //     //     cout << decoded_nh;
        //     //     s += 1;
        //     // }
        // }
        // cout << endl;

        return 0;
}

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
        if (mm > 1) return 0;
        if ((*bc_iter)[k] != bam_seqi(seq, j)) mm += 1;
        ++k;
    }
    return 1;
}
int compare_barcode(uint8_t* seq, vector<vector<int> >* barcodes, int start, int end, int* bc_idx) {
    vector<vector<int> >::iterator bc_iter = barcodes->begin();
    vector<vector<int> >::iterator bc_iter_end = barcodes->end();

    vector<int> mm(barcodes->size());

    int i = 0;
    int result = 0;
    for (; bc_iter != bc_iter_end ; ++bc_iter) {
        if (compare_barcode_local(bc_iter, seq, start, end))
            return i;
        ++i;
    }
    return -1;
    // // No perfect match found.
    // i = 0;
    // vector<int>::iterator mm_iter = mm.begin();
    // vector<int>::iterator mm_iter_end = mm.end();
    // int min_barcode = 1000;
    // for (; mm_iter != mm_iter_end; ++ mm_iter) {
    //     if (*mm_iter < min_barcode) {
    //         min_barcode = *mm_iter;
    //         *bc_idx = i;
    //     }
    //     ++i;
    // }
    // return min_barcode;
}

//intended for barcodes
int get_sequence(bam1_t* b, int start, int end, vector<vector<int> >* barcodes) {
    uint8_t* seq = bam_get_seq(b);
    uint8_t* qual = bam_get_qual(b);
    //uint32_t l_qseq = b->core.l_qseq;

    int min = 100000;
    int qual_int;
    for (int j = start; j <= end ; ++j) {
        int qual_int = int(*qual);
        if (qual_int < min) min = qual_int;
        qual += 1;
    }
    if (min < 20) return -1;

    int bc_idx = -1;
    int out = compare_barcode(seq, barcodes, start, end, &bc_idx);
    if (out > 1) return -1;
    return bc_idx;
}

int ShiftAdd(int sum, int digit)
{
     return sum*10 + digit;
}

//intended for umi
int get_sequence(bam1_t* b, int start, int end) {
    uint8_t* seq = bam_get_seq(b);
    uint8_t* qual = bam_get_qual(b);

    int min = 100000;
    int qual_int;
    for (int j = start; j <= end ; ++j) {
        int qual_int = int(*qual);
        if (qual_int < min) min = qual_int;
        qual += 1; //advance pointer
    }
    if (min < 20) {
        return 0;
    } else {
        vector<bitset<4> > umi_full(5);
        vector<int> umi_int(5);
        for (int j = start; j <= end; ++j) {
            umi_full[j] = bam_seqi(seq, j);
            umi_int[j] = umi_full[j].to_ulong();
        }
        return accumulate(umi_int.begin(), umi_int.end(), 0, ShiftAdd);
    }
}

int fill_db(BamDB* bamdb) {

    int tid = 0; // just first reference for now
    hts_itr_t* bam_itr = bam_itr_queryi(bamdb->get_idx(), tid, 0, bamdb->get_rlen(tid));
    dbRecord record = {};
    record.tid = tid;

    bam1_t* b = bam_init1();
    int result;
    vector<bitset<4> > umi_full(5);
    int i = 0;

    sqlite3_stmt * stmt;
    char * sErrMsg = 0;
    const char * tail = 0;
    char SQL[BUFFER_SIZE];


    sprintf(SQL, "INSERT INTO align VALUES (@IN, @FL, @CL, @TID, @POS, @STR, @BC, @UMI);");
    sqlite3_prepare_v2(bamdb->get_conn(),  SQL, BUFFER_SIZE, &stmt, &tail);


    sqlite3_exec(bamdb->get_conn(), "BEGIN TRANSACTION", NULL, NULL, &sErrMsg);

    while ((result = sam_itr_next(bamdb->get_bam(), bam_itr, b)) >= 0) {

        //if (i == 10) return 0;
        // get cigar
        if (bad_cigar(b)) continue;

        //const char tag[2] = {'N', 'H'};
        //filter_bad_reads(b, tag);

        split_qname(b, &record);

        if (bam_is_rev(b)) {
            record.strand = true;
            record.pos = bam_endpos(b);
        } else {
            record.pos = b->core.pos + 1;
        }

        record.bc = get_sequence(b, 5, 10, bamdb->get_barcodes());
        record.umi = get_sequence(b, 0, 4);

        try {
            insert_to_db(bamdb, &record, stmt);
        } catch ( exception& e ) {
            exit(1);
        }

        //record.reset_qname();
        ++i;
        // if (i == chunk_size) {
        //     commit
        // }
    }

    sqlite3_exec(bamdb->get_conn(), "END TRANSACTION", NULL, NULL, &sErrMsg);
    sqlite3_close(bamdb->get_conn());

    return 0;
}
