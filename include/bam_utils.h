#ifndef BAM_UTIL_H
#define BAM_UTIL_H

//#include <bam_sql_merge.h>
#include <string_utils.h>
#include <dbrecord.h>
#include <htslib/sam.h>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

uint64_t count_bam_records(hts_idx_t* idx, bam_hdr_t* header);

std::vector<int> seq2int(std::string& s);


int bad_cigar(bam1_t* b);

int filter_multi_reads(bam1_t* b);

int compare_barcode_local(vector<vector<int> >::iterator bc_iter, uint8_t* seq, int start, int end);

/*
Look for barcode match within given start and end of provided sequence.
Will also search sequence with position offsets provided in bc_offsets vector
*/
int compare_barcode(uint8_t* seq, vector<vector<int> >* barcodes, vector<int>* bc_offsets, \
                    int start, int end, int* used_offset);

/*
Find matching sequence in read[start:end] in barcodes vector
@barcodes, vector of int vectors, where nucleotides are encoded in 1,2,4,8 system
@bc_offsets, vector of positions offsets to try; will detect barcode if read indels
@min_qual, minimum quality score tolerated in barcode
@used_offset, variable set to bc_offset used for barcode ID
*/
int get_sequence(bam1_t* b, int start, int end, vector<vector<int> >* barcodes, \
                 vector<int>* bc_offsets, int min_qual, int* used_offset);
/*
Extract UMI sequence from read[start:end]
@used_offset, position offset used in barcode ID
*/
uint32_t get_sequence(bam1_t* b, int start, int end, int used_offset);

/*
used to combine vector of ints into one number
*/
int shift_add(int sum, int digit);

#endif
