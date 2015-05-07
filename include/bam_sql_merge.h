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
#include <map>
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
            //const char * dest_merge;
            // string dest_merge;
            // const char* dest_tmp;
            boost::filesystem::path dest_path;
            boost::filesystem::path merge_path;
            boost::filesystem::path tmp_path;
            const bool just_merge;
            const bool just_fill;
            samFile* bam;
            bam_hdr_t* header;
            map<int,char*> chroms;
            hts_idx_t* idx;
            map<int,sqlite3*> conns;
            sqlite3_stmt* stmt;

            vector<vector<int> > barcodes;
            vector<int> sequence_pos;
            vector<int> bc_offsets;

            int bc_min_qual;
    public:
            BamDB(const char* bam_fname, const char* final_fname, const char* barcodes_fname, int umi_length, int bc_min_qual, bool just_merge, bool just_fill);

            ~BamDB();

            /* Settors */
            void set_barcodes(const char* fname, vector<vector<int> >& vec_p);

            /* Gettors */

            samFile* get_bam() { return bam; }
            bam_hdr_t* get_header() { return header; }
            hts_idx_t* get_idx() { return idx; }
            int get_rlen(int tid) { return header->target_len[tid]; }
            char * get_chrom(int tid) { return chroms[tid]; }
            const boost::filesystem::path get_dest_path() { return dest_path; }
            const boost::filesystem::path get_tmp_path() { return tmp_path; }
            const boost::filesystem::path get_merge_path() { return merge_path; }

            sqlite3* get_conn() { return conns[0]; }
            sqlite3* get_conn(int index) { return conns[index]; }
            sqlite3_stmt* get_stmt() { return stmt; }

            int get_sequence_pos(int i) { return sequence_pos[i]; }
            vector<vector<int> >* get_barcodes() { return &barcodes; }
            vector<int>* get_bc_offsets() { return &bc_offsets; }
            int get_bc_min_qual() { return bc_min_qual; }
            bool merge_only() { return just_merge; }
            bool fill_only() { return just_fill; }

            /* Table creation */
            void create_align_table();
            void create_reftable();
            void create_rtree();
            void create_idcollapsed();

            /* Index creation */
            void index_cluster();
            void index_merge();
            void create_read1_index();
            void create_pos_indices();
            void create_collapsed_index();

            /* Processing */
            void merge_tables();
            void collapse_positions();
            void group_umi();

            /* Housekeeping */
            void remove_tmp_files();
            void drop_read_tables();
            void close_conn(int index);







    friend class dbRecord;
};

int create_index(BamDB* bamdb);
int fill_db(BamDB* bamdb);
/* template<typename TRecord> */
/* int fill_db_tid(BamDB* bamdb, TRecord* record, hts_itr_t* bam_itr); */

#endif //BAM_SQL_H
