#ifndef BAM_UTIL_H
#define BAM_UTIL_H

#include <string_utils.h>
//#include <dbrecord.h>
#include <htslib/sam.h>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <numeric>
#include <iostream>
using namespace std;

uint64_t count_bam_records(hts_idx_t* idx, bam_hdr_t* header);

std::vector<int> seq2int(std::string& s);

void bam_get_seq2(bam1_t* b, uint8_t* seq, uint8_t* qual);
void bam_get_seq_qual(bam1_t* b, uint8_t* seq, uint8_t* qual);
void bam_revcomp(bam1_t* b, int seqlen, uint8_t* seq_rc);
void bam_rev_qual(bam1_t* b, int seqlen, uint8_t* seq_rev);
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

int get_sequence(char* seq, char* qual, int start, int end, vector<const char*>* barcodes, \
                 vector<int>* bc_offsets, int min_qual, int* used_offset);
/*
Extract UMI sequence from read[start:end]
@used_offset, position offset used in barcode ID
*/
uint32_t get_sequence(bam1_t* b, int start, int end, int used_offset);
string get_sequence(char* seq, char* qual, int start, int end, int min_qual, int used_offset);

/*
used to combine vector of ints into one number
*/
int shift_add(int sum, int digit);
void print_uint8(uint8_t* arr, int seqlen, bool convert);

const char* convert_to_cstr(const string & s);

#endif
