#ifndef BAM_SQL_H
#define BAM_SQL_H

#include <vector>
#include <string.h>
#include <sqlite3.h>
#include <htslib/sam.h>

using namespace std;

#define BC_PATH "/home/brad/lib/barcodes/"
#define BUFFER_SIZE 256

//class BamDB;



class BamDB {
    private:
            const char* bam_fname;
            const char* dest_fname;
            const bool paired_end;
            samFile* bam;
            bam_hdr_t* header;
            hts_idx_t* idx;
            uint8_t pe;
            uint64_t total_mapped;
            uint64_t total_read;

            dbRecord* record;
            const char* db_fname;
            sqlite3* conn;
            sqlite3_stmt* stmt;

            vector<vector<int> > barcodes;
            vector<int> sequence_pos;
            vector<int> bc_offsets;

            int bc_min_qual;

    public:
            BamDB(const char* bam_fname, const char* dest_fname, const char* barcodes_fname, int umi_length, int bc_min_qual, bool paired_end);

            /* Settors */
            void set_barcodes(const char* fname, vector<vector<int> >& vec_p);

            /* Gettors */
            bool is_pe() { return paired_end; }
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
            sqlite3* get_conn() { return conn; }
            sqlite3_stmt* get_stmt() { return stmt; }

            int get_bc_min_qual() { return bc_min_qual; }
            dbRecord* get_record() { return record; }

            /* Other */
            void create_reftable();
            void increment_read() { ++total_read;}
            //void update_percent() { percent += 5; }
    friend class dbRecord;
};

class dbRecord : public BamDB {
    protected:
    //BamDB* bamdb;
    char instrument[BUFFER_SIZE];
    char flowcell[BUFFER_SIZE];
    char cluster[BUFFER_SIZE];
    int tid, bc, umi;
    vector<uint32_t> read_pos;
    bool strand;

public:
//dbRecord(BamDB* bamdb)
dbRecord()
{
//bamdb = bamdb;
};

/* Settors */
void set_instrument(const char* input) { strcpy(instrument, input); }
void set_flowcell(const char* input) { strcpy(flowcell, input); }
void set_cluster(const char* input) { strcpy(cluster, input); }
void set_tid(int tid) { tid = tid; }
void set_bc(BamDB* bamdb, bam1_t* b, int* used_offset);
void set_umi(bam1_t* b, int umi_start, int umi_end, int used_offset);

/* Gettors */
//const char* get_instrument() { return instrument; }
vector<uint32_t> get_pos() { return read_pos; }
virtual int set_positions(bam1_t* b, int read_num);
virtual int insert_to_db(sqlite3_stmt* stmt);
void split_qname(bam1_t* b);
};


class dbRecordSe: public dbRecord {
public:
//dbRecordSe(BamDB* bamdb)
dbRecordSe()
{
//dbRecord::bamdb = bamdb;
read_pos.resize(2);
};
int set_positions(bam1_t* b, int read_num) {
if (bam_is_rev(b)) {
read_pos[1] = b->core.pos + 1;
read_pos[0] = bam_endpos(b);
strand = true;
} else {
read_pos[0] = b->core.pos + 1;
read_pos[1] = bam_endpos(b);
strand = false;
}
return 0;
}
/* Other */
int insert_to_db(sqlite3_stmt* stmt);
};

class dbRecordPe: public dbRecord {
protected:
uint8_t insert;
public:
dbRecordPe() {
//dbRecord::bamdb = bamdb;
read_pos.resize(4);
};
// +1 offset for comparison with 1-based indexing
int set_positions(bam1_t* b, int read_num) {
if (read_num == 1) {
if (bam_is_rev(b)) {
read_pos[0] = bam_endpos(b);
read_pos[1] = b->core.pos + 1;
strand = true;
} else {
read_pos[0] = b->core.pos + 1;
read_pos[1] = bam_endpos(b);
strand = false;
}

} else if (read_num == 2) {
if (bam_is_rev(b)) {
read_pos[2] = bam_endpos(b);
read_pos[3] = b->core.pos + 1;
} else {
read_pos[2] = b->core.pos + 1;
read_pos[3] = bam_endpos(b);
}

}
return 0;
}

void set_insert(bam1_t* b) {
insert = b->core.isize;
}

/* Other */
int insert_to_db(sqlite3_stmt* stmt);
};



int create_table(BamDB* bamdb);
int create_index(BamDB* bamdb);
int fill_db(BamDB* bamdb);
int fill_db_tid(BamDB* bamdb, int tid, hts_itr_t* bam_itr);


#endif //BAM_SQL_H
