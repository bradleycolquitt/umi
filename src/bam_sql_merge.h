#ifndef BAM_SQL_H
#define BAM_SQL_H

#include <string_utils.h>
#include <bam_utils.h>
#include <sql_utils.h>
#include <dbrecord.h>
#include <sqlite_wrapper.h>
#include <timing.h>
#include <sqlite3.h>
#include <htslib/sam.h>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <libgen.h>
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
            char * dest_merge;
            char* dest_tmp;
            const bool paired_end;
            samFile* bam;
            bam_hdr_t* header;
            map<int,char*> chroms;
            hts_idx_t* idx;
            map<int,sqlite3*> conns;
            /* sqlite3* tmp_conn; */
            /* sqlite3* merge_conn; */
            sqlite3_stmt* stmt;

            vector<vector<int> > barcodes;
            vector<int> sequence_pos;
            vector<int> bc_offsets;

            int bc_min_qual;
    public:
            BamDB(const char* bam_fname, const char* final_fname, const char* barcodes_fname, int umi_length, int bc_min_qual, bool paired_end);
            ~BamDB() {
                sqlite3_close(conns[0]);
                sqlite3_close(conns[1]);
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
            char * get_chrom(int tid) { return chroms[tid]; }
            const char * get_dest_name() { return dest_fname; }
            char * get_tmp_name() { return dest_tmp; }
            sqlite3* get_conn() { return conns[0]; }
            sqlite3* get_conn(int index) { return conns[index]; }
            sqlite3_stmt* get_stmt() { return stmt; }

            int get_sequence_pos(int i) { return sequence_pos[i]; }
            vector<vector<int> >* get_barcodes() { return &barcodes; }
            vector<int>* get_bc_offsets() { return &bc_offsets; }
            int get_bc_min_qual() { return bc_min_qual; }

            /* Other */
                // Create 'reference table' containing human-readable reference names and tdd.
            void close_conn(int index);
            void create_reftable();
            int create_align_table();
            void index_cluster();
            int merge_tables();
            void drop_read_tables();
            void aggregate_umi();
            int create_rtree();
            int create_index();

    friend class dbRecord;
};

int create_index(BamDB* bamdb);
int fill_db(BamDB* bamdb);
/* template<typename TRecord> */
/* int fill_db_tid(BamDB* bamdb, TRecord* record, hts_itr_t* bam_itr); */

#endif //BAM_SQL_H
