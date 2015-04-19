#ifndef BAM_HASH_H
#define BAM_HASH_H

class BamRecord;

#include <bamrecord.h>
#include <string_utils.h>
#include <bam_utils.h>
#include <timing.h>
#include <htslib/sam.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>

using namespace std;

#define BC_PATH "/home/brad/lib/barcodes/"
#define BUFFER_SIZE 256

#ifdef DEBUG
#define DEBUG_LOG(x) cerr << x << endl;
#else
#define DEBUG_LOG(x)
#endif

class UmiHash
{
    public:
        void update(BamRecord * record);

        unordered_map<int, int >::iterator begin()
        {
            return umi_map.begin();
        }

        unordered_map<int, int >::iterator end()
        {
            return umi_map.end();
        }
    private:
        unordered_map<int, int> umi_map;
};

class BarcodeHash
{
    public:
        void update(BamRecord * record);

        unordered_map<int, unique_ptr<UmiHash> >::iterator begin()
        {
            return umi_map.begin();
        }

        unordered_map<int, unique_ptr<UmiHash> >::iterator end()
        {
            return umi_map.end();
        }

    private:
        unordered_map<int, unique_ptr<UmiHash> > umi_map;
};

class PositionHash
{
    public:
        void update(BamRecord * record);

        unordered_map<int, unique_ptr<BarcodeHash> >::iterator begin()
        {
            return barcode_map.begin();
        }

        unordered_map<int, unique_ptr<BarcodeHash> >::iterator end()
        {
            return barcode_map.end();
        }
    private:
        unordered_map<int, unique_ptr<BarcodeHash> > barcode_map;
};

class BamHash
{
    private:
            const char* bam_fname;
            const char* anno_fname;
            const char* dest_fname;
            boost::filesystem::path dest_path;

            samFile* bam;
            bam_hdr_t* header;
            map<int,char*> chroms;
            hts_idx_t* idx;

            vector<vector<int> > barcodes;
            vector<int> sequence_pos;
            vector<int> bc_offsets;
            int bc_min_qual;

            unordered_map<string, string> anno_map;
            unordered_map<string, unique_ptr<PositionHash> > position_map;

    public:
            BamHash(const char* bam_fname, const char* anno_name,
                    const char* final_fname, const char* barcodes_fname,
                    int umi_length, int bc_min_qual);

            ~BamHash();

            ofstream outfile;

            /* Settors */
            void set_barcodes(const char* fname, vector<vector<int> >& vec_p);

            /* Gettors */
            samFile* get_bam() { return bam; }
            bam_hdr_t* get_header() { return header; }
            hts_idx_t* get_idx() { return idx; }
            int get_rlen(int tid) { return header->target_len[tid]; }
            char * get_chrom(int tid) { return chroms[tid]; }
            const char * get_anno_fname() { return anno_fname; }
            const boost::filesystem::path get_dest_path() { return dest_path; }
            int get_sequence_pos(int i) { return sequence_pos[i]; }
            vector<vector<int> >* get_barcodes() { return &barcodes; }
            vector<int>* get_bc_offsets() { return &bc_offsets; }
            int get_bc_min_qual() { return bc_min_qual; }

            /* Other */
           void run();
           void hash_reads();
           void insert_anno(const string & read_id, const string & gene_id);
           bool find_read(bam1_t * b, BamRecord * record);
           void update_maps(BamRecord * record);
           void print_results();
};

void hash_annotation(BamHash* bamhash);
void hash_reads(BamHash* bamhash);
#endif //BAM_HASH_H
