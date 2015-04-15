#ifndef DBRECORD_H
#define DBRECORD_H

class BamDB;

#include <bam_sql_merge.h>
#include <string_utils.h>
#include <bam_utils.h>

#include <htslib/sam.h>
#include <sqlite3.h>
#include <iostream>
#include <vector>
#include <cstring>

#define BUFFER_SIZE 256

using namespace std;

class dbRecord0 {
    protected:
        BamDB* bamdb;
        char instrument[BUFFER_SIZE];
        char flowcell[BUFFER_SIZE];
        char cluster[BUFFER_SIZE];
        char chrom[BUFFER_SIZE];
        int32_t tid;
        vector<uint32_t> read_pos;
        bool strand;
    public:
        dbRecord0() {}
        dbRecord0(BamDB* bamdb)
            : bamdb(bamdb)
            , read_pos(2)
            , strand(false)
        {
        }
        virtual ~dbRecord0() {}

        /* Settors */
        int set_positions(bam1_t* b);
        void set_tid(int32_t tid);
        void set_chrom(BamDB* bamdb, int32_t tid);

        /* Gettors */
        const char* get_cluster() { return cluster; }

        /* Others */
        virtual int insert_to_db() = 0;
        void split_qname(bam1_t* b);

    friend class BamDB;
};

class dbRecord1: public dbRecord0 {
    private:
        int bc;
        uint32_t umi;
        //int32_t tid;
        sqlite3_stmt* insert_stmt;
    public:
        dbRecord1() {}
        dbRecord1(BamDB* bamdb);
        ~dbRecord1() {
            sqlite3_finalize(insert_stmt);
        }
        //void set_tid(int32_t tid) { tid = tid; }
        int32_t get_tid() { return tid; }
        void set_bc(BamDB* bamdb, bam1_t* b, int* used_offset);
        void set_umi(BamDB* bamdb, bam1_t* b, int used_offset);
        int insert_to_db();
};

class dbRecord2: public dbRecord0 {
    private:
        //int32_t tid;
        int32_t isize;
        sqlite3_stmt* insert_stmt;
    public:
        dbRecord2() {}
        dbRecord2(BamDB* bamdb);
        ~dbRecord2() {
            sqlite3_finalize(insert_stmt);
        }

        void set_isize(int32_t is);
        int insert_to_db();
};

#endif
