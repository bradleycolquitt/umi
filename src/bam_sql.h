#ifndef BAM_SQL_H
#define BAM_SQL_H

#include <vector>
#include <string>
#include <sqlite3.h>
#include <htslib/sam.h>

using namespace std;

#define BC_PATH "/home/brad/lib/barcodes/"
#define BUFFER_SIZE 256

class BamDB {
    private:
            const char* bam_fname;
            const char* dest_fname;
            samFile* bam;
            bam_hdr_t* header;
            hts_idx_t* idx;
            uint8_t pe;
            uint64_t total_mapped;
            uint64_t total_read;

            const char* db_fname;
            sqlite3* conn;
            vector<vector<int> > barcodes;
            vector<int> sequence_pos;
            vector<int> bc_offsets;

            int bc_min_qual;

    public:
            BamDB(const char* bam_fname, const char* dest_fname, const char* barcodes_fname, int umi_length, int bc_min_qual);
            void set_barcodes(const char* fname, vector<vector<int> >& vec_p);
            void create_reftable();
            void increment_read() { ++total_read;}
            //void update_percent() { percent += 5; }

            vector<int> get_sequence_pos() { return sequence_pos; }
            samFile* get_bam() { return bam; }
            bam_hdr_t* get_header() { return header; }
            hts_idx_t* get_idx() { return idx; }
            uint64_t get_mapped() { return total_mapped; }
            uint64_t get_read() { return total_read; }
            //int get_percent() { return percent_complete; }
            vector<vector<int> >* get_barcodes() { return &barcodes; }
            vector<int>* get_bc_offsets() { return &bc_offsets; }
            int get_rlen(int tid) { return header->target_len[tid]; }
            sqlite3* get_conn() {return conn;}
            int get_bc_min_qual() { return bc_min_qual; }

};

struct dbRecord {
    char instrument[BUFFER_SIZE];
    char flowcell[BUFFER_SIZE];
    char cluster[BUFFER_SIZE];

    int tid = 0;
    uint32_t pos_head = 0;
    uint32_t pos_tail;
    bool strand;
    int bc;
    int umi;

};

int create_table(BamDB* bamdb);
int create_index(BamDB* bamdb);
int fill_db(BamDB* bamdb);
int fill_db_tid(BamDB* bamdb, int tid, hts_itr_t* bam_itr);


#endif //BAM_SQL_H
