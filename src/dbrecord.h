#ifndef DBRECORD_H
#define DBRECORD_H

class BamDB;

#include <bam_sql_serial2.h>
#include <string_utils.h>
#include <bam_utils.h>

#include <htslib/sam.h>
#include <sqlite3.h>
#include <iostream>
#include <vector>
#include <cstring>

#define BUFFER_SIZE 256

using namespace std;

class dbRecord {
    protected:
        BamDB* bamdb;
        char instrument[BUFFER_SIZE];
        char flowcell[BUFFER_SIZE];
        char cluster[BUFFER_SIZE];
        int tid, bc, umi;
        vector<uint32_t> read_pos;
        bool strand;
        sqlite3_stmt* stmt_exists;

    public:
        //dbRecord() {}
        dbRecord(BamDB* bamdb);
        /* dbRecord(BamDB* bamdb) */
        /*     : bamdb(bamdb) */
        /* { */
        /*    prepare_exists(); */
        /* } */
        virtual ~dbRecord();

        /* Settors */
        void set_tid(int tid) { tid = tid; }
        void set_bc(BamDB* bamdb, bam1_t* b, int* used_offset);
        void set_umi(BamDB* bamdb, bam1_t* b, int used_offset);
        int set_positions(bam1_t* b);
        virtual int set_positions(bam1_t* b, int read_num) = 0;

        /* Gettors */
        const char* get_cluster() { return cluster; }

        /* Others */
        void prepare_exists();
        virtual int insert_to_db() = 0;

        void split_qname(bam1_t* b);

    friend class BamDB;
};

class dbRecordSe: public dbRecord {
    protected:
        sqlite3_stmt* stmt;
    public:
        dbRecordSe(BamDB* bamdb);
        int set_positions(bam1_t* b, int read_num);
        /* Other */
        int insert_to_db();
        friend class BamDB;
};

class dbRecordPe: public dbRecord {
    protected:
        sqlite3_stmt* insert_stmt;
        sqlite3_stmt* update_stmt1;
        sqlite3_stmt* update_stmt2;
        uint32_t insert;
    public:
        dbRecordPe(BamDB* bamdb);

        // +1 offset for comparison with 1-based indexing
        int set_positions(bam1_t* b, int read_num);
        void set_insert(bam1_t* b) { insert = abs(b->core.isize); }

        /* Other */
        int exists(BamDB* bamdb);
        int insert_to_db();
        int insert_to_db(int read_num);
        int update_record(int read_num);
    friend class BamDB;
};

class dbRecord0 {
    protected:
        BamDB* bamdb;
        char instrument[BUFFER_SIZE];
        char flowcell[BUFFER_SIZE];
        char cluster[BUFFER_SIZE];
        vector<uint32_t> read_pos;
        bool strand;
    public:
        dbRecord0() {}
        dbRecord0(BamDB* bamdb)
            : bamdb(bamdb)
            , read_pos(2)
        {
        }
        virtual ~dbRecord0() {}

        /* Settors */
        int set_positions(bam1_t* b);

        /* Gettors */
        const char* get_cluster() { return cluster; }

        /* Others */
        virtual int insert_to_db() = 0;
        void split_qname(bam1_t* b);

    friend class BamDB;
};

class dbRecord1: public dbRecord0 {
    private:
        int bc, umi;
        int32_t tid;
        sqlite3_stmt* insert_stmt;
    public:
        dbRecord1() {}
        dbRecord1(BamDB* bamdb);
        void set_tid(int32_t tid) { tid = tid; }
        int32_t get_tid() { return tid; }
        void set_bc(BamDB* bamdb, bam1_t* b, int* used_offset);
        void set_umi(BamDB* bamdb, bam1_t* b, int used_offset);
        int insert_to_db();
};

class dbRecord2: public dbRecord0 {
    private:
        int32_t tid;
        sqlite3_stmt* insert_stmt;
    public:
        dbRecord2() {}
        dbRecord2(BamDB* bamdb);
        void set_tid(int32_t tid) { tid = tid; }
        int insert_to_db();
};

#endif
