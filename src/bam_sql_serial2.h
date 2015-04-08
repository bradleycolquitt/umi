#ifndef BAM_SQL_H
#define BAM_SQL_H

#include <string_utils.h>
#include <bam_utils.h>
#include <dbrecord.h>
#include <sqlite3.h>
#include <htslib/sam.h>
#include <boost/regex.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

#define BC_PATH "/home/brad/lib/barcodes/"
#define BUFFER_SIZE 256

#ifdef DEBUG
#define DEBUG_LOG(x) cerr << x << endl;
#else
#define DEBUG_LOG(x)
#endif


/*
Iter through Bam
    Insert into read1/2 dbs
Join read1/2 dbs on cluster
    Insert into one db
*/

class BamDB {
    private:
            const char* bam_fname;
            const char* dest_fname;
            //vector<const char*> dest_fname;
            /* const char* dest_fname1; //stores read1 */
            /* const char* dest_fname2; //stores read2 */
            /* const char* final_fname; //final merged db */
            const bool paired_end;
            samFile* bam;
            bam_hdr_t* header;
            hts_idx_t* idx;
            //vector<sqlite3*> conns;
            /* sqlite3* conn1;          //connection to read1 db */
            /* sqlite3* conn2;          //connection to read2 db */
            sqlite3* conn;           //connection to merge db
            sqlite3_stmt* stmt;

            vector<vector<int> > barcodes;
            vector<int> sequence_pos;
            vector<int> bc_offsets;

            int bc_min_qual;
    public:
            BamDB(const char* bam_fname, const char* final_fname, const char* barcodes_fname, int umi_length, int bc_min_qual, bool paired_end);
            ~BamDB() {
//                close_connections();
                sqlite3_close(conn);
                sam_close(bam);
                bam_hdr_destroy(header);
                hts_idx_destroy(idx);
            }

            /* Settors */
            void set_barcodes(const char* fname, vector<vector<int> >& vec_p);

            /* Gettors */
            bool is_pe() { return paired_end; }
            samFile* get_bam() { return bam; }
            bam_hdr_t* get_header() { return header; }
            hts_idx_t* get_idx() { return idx; }
            int get_rlen(int tid) { return header->target_len[tid]; }
            sqlite3* get_conn() { return conn; }
            sqlite3_stmt* get_stmt() { return stmt; }

            int get_sequence_pos(int i) { return sequence_pos[i]; }
            vector<vector<int> >* get_barcodes() { return &barcodes; }
            vector<int>* get_bc_offsets() { return &bc_offsets; }
            int get_bc_min_qual() { return bc_min_qual; }

            /* Other */
            // Create 'reference table' containing human-readable reference names and tid.
            void create_reftable();
            int create_align_table();
            void index_cluster();

            /* void close_connections() { */
            /*     for (vector<sqlite*>::iterator iter=conns.begin(); iter != conns.end; ++iter) { */
            /*         sqlite3_close(*iter); */
            /*     } */
            /* } */
    friend class dbRecord;
};



//int create_table(BamDB* bamdb);
int create_index(BamDB* bamdb);
int fill_db(BamDB* bamdb);
template<typename TRecord>
int fill_db_tid(BamDB* bamdb, TRecord* record, hts_itr_t* bam_itr);

#endif //BAM_SQL_H
