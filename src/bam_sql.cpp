#include <fstream>
#include <bam_sql.h>
#include <bam_utils.h>
#include <string_utils.h>
#include <boost/regex.hpp>

//#define DEBUG

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

static int callback_table(void *NotUsed, int argc, char **argv, char **azColName){
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
    const char* statement = "CREATE TABLE IF NOT EXISTS align (name text,       \
                                                           instrument text,     \
                                                           flowcell text,       \
                                                           chrom int,           \
                                                           position int,        \
                                                           strand int,          \
                                                           bc int,              \
                                                           umi text);";

    rc = sqlite3_exec(bamdb->get_conn(), statement, callback_table, 0, &err_msg);
    if( rc != SQLITE_OK ){
       fprintf(stderr, "SQL error: %s\n", err_msg);
       sqlite3_free(err_msg);
    } else{
       fprintf(stdout, "Align table created successfully\n");
    }
    return 0;
}

int process_cigar(bam1_t* b, uint32_t* cigar) {
    cigar = bam_get_cigar(b);
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

void split_qname(bam1_t* b, char** instrument, char** flowcell) {
    char* qname = bam_get_qname(b);
    char* qname_array = strtok(qname, ":");
    int j = 0;
    while (qname_array != NULL) {
        if (j == 0) {
            *instrument = qname_array;
        } else if (j == 2) {
            *flowcell = qname_array;
        }
        qname_array = strtok(NULL, ":");
        ++j;
    }
}

int compare_barcode(uint8_t* seq, vector<vector<int> >* barcodes, int start, int end, int* bc_idx) {
    vector<vector<int> >::iterator bc_iter = barcodes->begin();
    vector<vector<int> >::iterator bc_iter_end = barcodes->end();

    vector<int> mm(barcodes->size());

    int i = 0;
    int k = 0;
    for (; bc_iter != bc_iter_end ; ++bc_iter) {
        for (int j = start; j <= end ; ++j) {
            #ifdef DEBUG
            cout << (*bc_iter)[k] << " ";
            cout << "i" << bam_seqi(seq, j) << " ";
            #endif
            if ((*bc_iter)[k] != bam_seqi(seq, j)) mm[k] += 1;
            ++k;
        }
        k = 0;
        #ifdef DEBUG
        cout << endl;
        #endif
        if (mm[i] == 0) {
            *bc_idx = i;
            return 0;
        }
        ++i;
    }

    // No perfect match found.
    i = 0;
    vector<int>::iterator mm_iter = mm.begin();
    vector<int>::iterator mm_iter_end = mm.end();
    int min_barcode = 1000;
    for (; mm_iter != mm_iter_end; ++ mm_iter) {
        if (*mm_iter < min_barcode) {
            min_barcode = *mm_iter;
            *bc_idx = i;
        }
        ++i;
    }
    return min_barcode;
}
int get_barcode(bam1_t* b, int start, int end, vector<vector<int> >* barcodes) {
    uint8_t* seq = bam_get_seq(b);
    uint8_t* qual = bam_get_qual(b);
    uint32_t l_qseq = b->core.l_qseq;

    int min = 100000;
    int qual_int;
    for (int j = 0; j < b->core.l_qseq ; ++j) {
        int qual_int = int(*qual);
        if (qual_int < min) min = qual_int;
        qual += 1;
    }
    if (min < 20) return -1;

    int bc_idx = -1;
    int out = compare_barcode(seq, barcodes, start, end, &bc_idx);
    //cout << bc_idx << endl;



}

int fill_db(BamDB* bamdb) {
    //boost::regex filter_pattern ( "[IDN]" );
    int tid = 0; // just first reference for now
    hts_itr_t* bam_itr = bam_itr_queryi(bamdb->get_idx(), tid, 0, bamdb->get_rlen(tid));

    uint32_t* cigar;
    uint8_t* seq;
    uint8_t* qual;
    uint8_t* aux_p;

    char* instrument;
    char* flowcell;

    uint32_t pos;
    // while (bam_itr_next(bamdb->get_bam(), bam_itr, r)) {
    // }


    bam1_t* b = bam_init1();
    int result;
    int i = 0;
    while ((result = sam_itr_next(bamdb->get_bam(), bam_itr, b)) >= 0) {
        if (i == 10) return 0;
        // get cigar
        if (process_cigar(b, cigar)) {
            continue;
        }


        const char tag[2] = {'N', 'H'};
        //filter_bad_reads(b, tag);


        if (i==0) { split_qname(b, &instrument, &flowcell); }

        if (bam_is_rev(b)) {
            pos = bam_endpos(b);
        } else {
            pos = b->core.pos + 1;
        }

        seq = bam_get_seq(b);
        qual = bam_get_qual(b);

        for (int j = 0; j < b->core.l_qseq ; ++j) {
            //cout << int(*qual) << " ";
            qual += 1;
        }
        //cout << endl;
        get_barcode(b, 5, 10, bamdb->get_barcodes());
        ++i;

    }




//         qname = *bam_get_qname(r);
//         vector<string> qname_list = split(qname); // need to include from R project
//         instrument = qname_list[0];
//         flowcell = qname_list[2];

//         if (bam_is_rev(r)) {
//             pos = bam_endpos(r);
//         } else {
//             pos = bam->core.pos;
//         }

//         seq_p = bam_get_seq(r);
//         qual_p = bam_get_qual(r);

//         bc = get_barcode(seq_p, qual_p, 5, 11, bamdb->barcodes);
//         umi = get_umi(seq_p, qual_p, 0, 4);
//     }

     return 0;
}


void get_umi(uint32_t* seq_p, uint32_t* qual_p, int start, int end, string& umi) {

}
