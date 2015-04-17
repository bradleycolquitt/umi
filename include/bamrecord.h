#ifndef BAMRECORD_H
#define BAMRECORD_H

class BamHash;

#include <bam_hash.h>
#include <htslib/sam.h>
#include <string>

using namespace std;

class BamRecord
{
    private:
        string read_id;
        string gene_id;
        int32_t tid;
        uint32_t pos;
        int bc;
        uint32_t umi;
    public:
        BamRecord(int tid);
        string get_read_id() { return read_id; }
        string get_gene_id() { return gene_id; }
        uint32_t get_pos() { return pos; }
        int get_bc() { return bc; }
        uint32_t get_umi() { return umi; }

        void set_read_id( string val ) { read_id = val; }
        void set_gene_id( string val ) { gene_id = val; }
        void set_position( uint32_t val ) { pos = val; }
        void set_bc(BamHash * bamhash, bam1_t * b, int * used_offset);
        void set_umi(BamHash * bamhash, bam1_t * b, int used_offset);
};

#endif
