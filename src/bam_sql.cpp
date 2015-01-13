#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <stdexcept>
#include <sqlite3.h>
#include <bam_sql.h>
#include <htslib/sam.h>
#include <bam_utils.h>
#include <boost/regex.hpp>
using namespace std;

#define BC_PATH "/home/brad/lib/barcodes/"

// class BamDB {

//     // Variables
//     private:
//         const char* bam_fname;
//         const char* dest_fname;
//         samFile* bam;
//         bam_hdr_t* header;
//         hts_idx_t* idx;
//         uint64_t total_mapped;
//         const char* db_fname;
//         sqlite3* conn;
//         vector<string> barcodes;

//     public:
//         // Settors
//         void set_barcodes(const char* fname, vector<string>& vec_p);
//         // Gettors
//         hts_idx_t* get_idx() { return idx; }
//         vector<string>* get_barcodes() { return &barcodes; }
//         int get_rlen(int tid) { return header->target_len[tid]; }
// };

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
        cout << bc_path << endl;
        try {
            set_barcodes(bc_path, barcodes);
        } catch (exception &e) {
            cout << e.what() << endl;
        }
    }

// Read in barcodes file and load sequences intobarcodes vector
// consider encoding barcodes as bits or just int
void BamDB::set_barcodes(const char* fname, vector<string>& vec_p) {

    ifstream bc_s(fname, ifstream::in);

    if (!bc_s.good()) {
         throw runtime_error("Error: Barcode file not found.");
    }

    string line;
    while(getline(bc_s, line)) {
        vec_p.push_back(line);
    }
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

int fill_db(BamDB* bamdb) {
    boost::regex pattern ( "[IDN]" );
    int tid = 0; // just first reference for now
    hts_itr_t* bam_itr = bam_itr_queryi(bamdb->get_idx(), tid, 0, bamdb->get_rlen(tid));
    bam1_t* r;

    uint32_t* cigar;
    uint8_t* aux_p;
    char* aux;

    string qname;
    int use = 1;

    uint32_t pos;
    while (bam_itr_next(bamdb->get_bam(), bam_itr, r)) {
    }

//         // get cigar
//         cigar = bam_get_cigar(r);

//         for (int i=0; i<r->core.n_cigar; ++i) {
//             uint32_t op = bam_cigar_op(cigar[i]);
//             if (op > 0 & op < 4) {
//                 use = 0;
//                 break;
//             }
//         }
//         if (!use) continue;

//         // get aux. filter out NH>1
//         aux_p = bam_get_aux(r);
//         aux = bam_aux2A(aux_p);
//         size_t aux_len = strlen(aux);

//         boost::regex pattern = "(NH\\:[0-9]+)";
//         boost::match_results<string::const_iterator> what;
//         while(boost::regex_search(aux, what, pattern, flags)) {
//             char* match = what[0];
//             if (strlen(match) > 4) {
//                 use = 0;
//                 break;
//             }
//             int mm = atoi(match[3]);
//             if (mm > 1) {
//                 use = 0;
//                 break;
//             }
//         }
//         if (!use) continue;

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

int get_barcode(uint32_t* seq_p, uint32_t* qual_p, int start, int end, vector<string>& barcodes) {}
void get_umi(uint32_t* seq_p, uint32_t* qual_p, int start, int end, string& umi) {

}
