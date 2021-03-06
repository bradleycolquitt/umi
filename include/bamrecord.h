#ifndef BAMRECORD_H
#define BAMRECORD_H

class BamHash;

#include <bam_hash2.h>
#include <htslib/sam.h>
#include <string>
#include <iostream>
#include <memory>
using namespace std;

class BamRecord
{
    private:
        string read_id;
        string gene_id;
        int32_t tid;
        bool strand;
        uint32_t pos;
        //std::shared_ptr<bam1_t> brs;
        bam1_t* br;
        int bc;
        uint32_t umi;
        string umi2;
    public:
        BamRecord(int tid);
        BamRecord(int tid, int pos, string read_id);
        BamRecord(bam1_t* br, string read_id);
        BamRecord(const string read_id, const string gene_id);
        string get_read_id() { return read_id; }
        string get_gene_id() { return gene_id; }
        bool get_strand() { return strand; }
        int32_t get_tid() { return tid; }
        uint32_t get_pos() { return pos; }
        int get_bc() { return bc; }
        uint32_t get_umi() { return umi; }
        string get_umi2() { return umi2; }
        bam1_t* get_bamread() { return br; }

        void set_read_id( string val ) { read_id = val; }
        void set_gene_id( string val ) { gene_id = val; }
        void set_strand ( bam1_t * b) { strand = (b->core.flag&BAM_FREVERSE) == 0; }
        void set_tid(uint32_t val) { tid = val; }
        void set_position( uint32_t val ) { pos = val; }
        //void set_br( shared_ptr<bam1_t> b) { br = b; }
        void set_bc(BamHash * bamhash, bam1_t * b, int * used_offset);
        void set_umi(BamHash * bamhash, bam1_t * b, int used_offset);
        void set_bc(BamHash * bamhash, char * seq, char * qual, int * used_offset);
        void set_umi(BamHash * bamhash, char * seq, char * qual, int used_offset);

        //void set_bam1(BamHash * bamhash);
        bool is_complete();
};

#endif
