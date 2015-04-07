#ifndef BAM_SQL_H
#define BAM_SQL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <sqlite3.h>
#include <htslib/sam.h>

using namespace std;

#define BC_PATH "/home/brad/lib/barcodes/"
#define BUFFER_SIZE 256

class dbRecord;

class BamDB {
  private:
            const char* bam_fname;
            const char* dest_fname;
            const bool paired_end;
            samFile* bam;
            bam_hdr_t* header;
            hts_idx_t* idx;
            sqlite3* conn;
            sqlite3_stmt* stmt;

            vector<vector<int> > barcodes;
            vector<int> sequence_pos;
            vector<int> bc_offsets;

            int bc_min_qual;

    public:
            BamDB(const char* bam_fname, const char* dest_fname, const char* barcodes_fname, int umi_length, int bc_min_qual, bool paired_end);
            ~BamDB() {
                cerr << "destroy BamDB" << endl;
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
            void create_reftable();
            int create_align_table();

    friend class dbRecord;
};

class dbRecord {
    protected:
        BamDB* bamdb;
        char instrument[BUFFER_SIZE];
        char flowcell[BUFFER_SIZE];
        char cluster[BUFFER_SIZE];
        int tid, bc, umi;
        vector<uint32_t> read_pos;
        bool strand;

    public:
        dbRecord() {}
        dbRecord(BamDB* bamdb);
        virtual ~dbRecord() {}
        /* Settors */
        /* void set_instrument(const char* input) { strcpy(instrument, input); } */
        /* void set_flowcell(const char* input) { strcpy(flowcell, input); } */
        /* void set_cluster(const char* input) { strcpy(cluster, input); } */
        void set_tid(int tid) { tid = tid; }
        void set_bc(BamDB* bamdb, bam1_t* b, int* used_offset);
        void set_umi(BamDB* bamdb, bam1_t* b, int used_offset);
        virtual int set_positions(bam1_t* b, int read_num) = 0;

        /* Gettors */
        const char* get_cluster() { return cluster; }
        vector<uint32_t> get_pos() { return read_pos; }

        virtual vector<int> insert_to_db() = 0;
        void split_qname(bam1_t* b);

    friend class BamDB;
};

class dbRecordSe: public dbRecord {
    protected:
        sqlite3_stmt* stmt;
    public:
        dbRecordSe(BamDB* bamdb) {
            read_pos.resize(2);

            //Prepare insert statement
            const char* tail = 0;
            char SQL[BUFFER_SIZE];
            sprintf(SQL, "INSERT INTO align VALUES (@IN, @FL, @CL, @TID, @HPOS, @TPOS, @STR, @BC, @UMI);");
            sqlite3_prepare_v2(bamdb->get_conn(), SQL, BUFFER_SIZE, &stmt, &tail);

        }

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
        vector<int> insert_to_db();
        friend class BamDB;
};

class dbRecordPe: public dbRecord {
    protected:
        sqlite3_stmt* insert_stmt;
        sqlite3_stmt* insert_stmt1;
        sqlite3_stmt* insert_stmt2;
        sqlite3_stmt* update_stmt1;
        sqlite3_stmt* update_stmt2;
        uint8_t insert;
    public:
        dbRecordPe(BamDB* bamdb) {
            read_pos.resize(4);
            const char* tail = 0;

            char insert_sql[BUFFER_SIZE];
            sprintf(insert_sql, "INSERT INTO align VALUES (@IN, @FL, @CL, @TID, @HPOS1, @TPOS1, @HPOS2, @TPOS2, @INS, @STR, @BC, @UMI);");
            cerr << "prep insert1" << endl;
            sqlite3_prepare_v2(bamdb->get_conn(), insert_sql, BUFFER_SIZE, &insert_stmt, &tail);
            //Prepare insert1 statement
            char insert_sql1[BUFFER_SIZE];
            sprintf(insert_sql1, "INSERT INTO align VALUES (@IN, @FL, @CL, @TID, @HPOS1, @TPOS1, @STR, @BC, @UMI);");
            cerr << "prep insert1" << endl;
            sqlite3_prepare_v2(bamdb->get_conn(), insert_sql1, BUFFER_SIZE, &insert_stmt1, &tail);

            //Prepare insert2 statement
            char insert_sql2[BUFFER_SIZE];
            sprintf(insert_sql2, "INSERT INTO align VALUES (@HPOS2, @HPOS2, @INS);");
            cerr << "prep insert2" << endl;
            sqlite3_prepare_v2(bamdb->get_conn(), insert_sql2, BUFFER_SIZE, &insert_stmt2, &tail);

            //Prepare update statement, read1
            char update_sql1[BUFFER_SIZE];
            sprintf(update_sql1, "UPDATE align SET @HPOS1=%d, @TPOS1=%d, @STR=%d, @BC=%d, @UMI=%d");
            cerr << "prep update1" << endl;
            sqlite3_prepare_v2(bamdb->get_conn(), update_sql1, BUFFER_SIZE, &update_stmt1, &tail);
            //Prepare update statement, read2
            char update_sql2[BUFFER_SIZE];
            sprintf(update_sql2, "UPDATE align SET @HPOS2=%d, @TPOS2=%d, @INS=%d");
            cerr << "prep update2" << endl;
            sqlite3_prepare_v2(bamdb->get_conn(), update_sql2, BUFFER_SIZE, &update_stmt2, &tail);
        }

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
        void set_insert(bam1_t* b) { insert = b->core.isize; }

        /* Other */
        vector<int> insert_to_db();
        vector<int> insert_to_db(int read_num);
        int update_record(int read_num);
    friend class BamDB;
};

int create_table(BamDB* bamdb);
int create_index(BamDB* bamdb);
int fill_db(BamDB* bamdb);
template<typename TRecord>
int fill_db_tid(BamDB* bamdb, TRecord* record, hts_itr_t* bam_itr);

#endif //BAM_SQL_H
