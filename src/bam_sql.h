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
            uint64_t total_mapped;
            const char* db_fname;
            sqlite3* conn;
            vector<vector<int> > barcodes;

    public:
            BamDB(const char* bam_fname, const char* dest_fname, const char* barcodes_fname);
            void set_barcodes(const char* fname, vector<vector<int> >& vec_p);

            samFile* get_bam() { return bam; }
            bam_hdr_t* get_header() { return header; }
            hts_idx_t* get_idx() { return idx; }
            vector<vector<int> >* get_barcodes() { return &barcodes; }
            int get_rlen(int tid) { return header->target_len[tid]; }
            sqlite3* get_conn() {return conn;}
};

struct dbRecord {
    char instrument[BUFFER_SIZE];
    char flowcell[BUFFER_SIZE];
    char cluster[BUFFER_SIZE];

    int tid;
    uint32_t pos;
    bool strand;
    int bc;
    int umi;

};

int create_table(BamDB* bamdb);

int fill_db(BamDB* bamdb);

#endif //BAM_SQL_H
