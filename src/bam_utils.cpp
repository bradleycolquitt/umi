#include <bam_utils.h>

using namespace std;

uint64_t count_bam_records(hts_idx_t* idx, bam_hdr_t* header) {

    uint64_t* mapped = new uint64_t[header->n_targets * sizeof(uint64_t)];
    uint64_t* unmapped = new uint64_t[header->n_targets * sizeof(uint64_t)];

    int* ret;
    ret = new int[header->n_targets];

    for (int i=0; i < header->n_targets; ++i) {
        ret[i] = hts_idx_get_stat(idx, i, &mapped[i], &unmapped[i]);
    }

    int total = 0;
    for (int i=0; i < header->n_targets; ++i) {
        total += mapped[i];
    }
    return total;
}

vector<int> seq2int(string& s) {
    string::iterator it = s.begin();
    string::iterator it_end = s.end();
    vector<int> out_vec(s.size());
    for (int i=0; i < s.size() ; ++i) {

    switch ( s[i] )
    {
    case 'A':
        out_vec[i] = 1;
        break;
    case 'C':
        out_vec[i] = 2;
        break;
    case 'G':
        out_vec[i] = 4;
        break;
    case 'T':
        out_vec[i] = 8;
        break;
    case 'N':
        out_vec[i] = 15;
        break;
    }
    }
    return out_vec;
}

void bam_get_seq2(bam1_t* b, int seqlen, uint8_t* seq_p, uint8_t* qual_p) {
    uint8_t* seq = bam_get_seq(b);
    for (int i = 0; i < seqlen; ++i) {
        seq_p[i] = bam_seqi(seq, i);
    }
    memcpy(qual_p, bam_get_qual(b), seqlen);
}

void bam_get_seq_qual(bam1_t* b, uint8_t* seq_p, uint8_t* qual_p) {
    int seqlen = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
    if (bam_is_rev(b)) {
        //bam_get_seq2(b, seqlen, seq_p, qual_p);
        bam_revcomp(b, seqlen, seq_p);
        bam_rev_qual(b, seqlen, qual_p);
    } else  {
        bam_get_seq2(b, seqlen, seq_p, qual_p);
    }
}

void bam_revcomp(bam1_t* b, int seqlen, uint8_t* seq_rc) {
    uint8_t* seq = bam_get_seq(b);
    for (int i = 0 ; i < seqlen ; ++i) {
        switch(bam_seqi(seq, i)) {
           case 1 :
               seq_rc[seqlen - i - 1] = 8;
               break;
           case 2 :
               seq_rc[seqlen - i - 1] = 4;
               break;
           case 4 :
               seq_rc[seqlen - i - 1] = 2;
               break;
           case 8 :
               seq_rc[seqlen - i - 1] = 1;
               break;
           case 15 :
               seq_rc[seqlen - i - 1] = 15;
               break;
        }
    }
}

void bam_rev_qual(bam1_t* b, int seqlen, uint8_t* seq_rev) {
    uint8_t* qual = bam_get_qual(b);
    for (int i = 0 ; i < seqlen ; ++i) {
        seq_rev[seqlen - i - 1] = int(qual[i]);
    }
}

/***********************
Read Processing utilities
**************************/
int bad_cigar(bam1_t* b) {
    uint32_t* cigar = bam_get_cigar(b);
    int op;
    for (int i=0; i < b->core.n_cigar; ++i) {
        op = bam_cigar_op(cigar[i]);
        if (op == BAM_CINS || op == BAM_CDEL || op == BAM_CREF_SKIP) {
            return 1;
        }
    }
    return 0;
}

// Filter out multimappers (NH>1)
int filter_multi_reads(bam1_t* b) {
    uint8_t *s = bam_aux_get(b, "NH");
    if (s) {
        if (bam_aux2i(s) > 1) {
            return 1;
        }
    }
    return 0;
}

int compare_barcode_local(vector<vector<int> >::iterator bc_iter, uint8_t* seq, int start, int end){
    int k = 0;
    int mm = 0;

    for (int j = start; j <= end ; ++j) {
        if (mm > 1) return 0; // return after second mismatch
        if ((*bc_iter)[k] != seq[j]) mm += 1;
        ++k;
    }
    return 1; // match found
}

int compare_barcode_local(vector<const char*>::iterator bc_iter, char* seq, int start, int end)
{
    int k = 0;
    int mm = 0;

    for (int j = start; j <= end ; ++j)
    {
        if (mm > 1) return 0; // return after second mismatch

        if ((*bc_iter)[k] != seq[j]) mm += 1;
        ++k;
    }

    return 1; // match found
}

/*
Look for barcode match within given start and end of provided sequence.
Will also search sequence with position offsets provided in bc_offsets vector
*/
int compare_barcode(uint8_t* seq, vector<vector<int> >* barcodes, vector<int>* bc_offsets, int start, int end, int* used_offset) {

    vector<int>::iterator offsets_iter = bc_offsets->begin();
    vector<int>::iterator offsets_iter_end = bc_offsets->end();
    vector<vector<int> >::iterator bc_iter;
    vector<vector<int> >::iterator bc_iter_end = barcodes->end();

    int i = 0;
    for (; offsets_iter != offsets_iter_end; ++offsets_iter) {
        bc_iter = barcodes->begin();
        for (; bc_iter != bc_iter_end ; ++bc_iter) {

            if (compare_barcode_local(bc_iter, seq, start + *offsets_iter, end + *offsets_iter)) {

                *used_offset = *offsets_iter;
                return i;
            }
            ++i;
        }
        i = 0;
    }
    return -2;
}

int compare_barcode(char* seq, vector<const char*>* barcodes, vector<int>* bc_offsets, int start, int end, int* used_offset) {

    vector<int>::iterator offsets_iter = bc_offsets->begin();
    vector<int>::iterator offsets_iter_end = bc_offsets->end();
    vector<const char*>::iterator bc_iter;
    vector<const char*>::iterator bc_iter_end = barcodes->end();

    int i = 0;
    for (; offsets_iter != offsets_iter_end; ++offsets_iter) {
        bc_iter = barcodes->begin();
        for (; bc_iter != bc_iter_end ; ++bc_iter) {

            if (compare_barcode_local(bc_iter, seq, start + *offsets_iter, end + *offsets_iter)) {
                *used_offset = *offsets_iter;
                return i;
            }
            ++i;
        }
        i = 0;
    }
    return -2;
}

//for barcodes
int get_sequence(bam1_t* b, int start, int end, vector<vector<int> >* barcodes, \
                 vector<int>* bc_offsets, int min_qual, int* used_offset)
{
    int seqlen = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
    uint8_t seq[seqlen];
    uint8_t qual[seqlen];
    bam_get_seq_qual(b, seq, qual);

    // skip read if barcode quality less than min_qual
    int min = 100000;
    int qual_int;
    for (int j = start; j <= end ; ++j) {
        qual_int = int(qual[j]);
        if (qual_int < min) min = qual_int;
    }
    if (min < min_qual) return -1;

    // returns index of perfect match or one mismatch
    return compare_barcode(seq, barcodes, bc_offsets, start, end, used_offset);
}

int get_sequence(char* seq, char* qual, int start, int end, vector<const char*>* barcodes, \
                 vector<int>* bc_offsets, int min_qual, int* used_offset)
{
    //cout << "here get_sequence" << endl;

    int min = 100000;
    int qual_int;
    for (int j = start; j <= end ; ++j) {
        qual_int = int(qual[j]) - 33;
        //cout << qual[j] << " " << qual_int << endl;
        if (qual_int < min) min = qual_int;
    }
    if (min < min_qual) return -1;

    // returns index of perfect match or one mismatch; +1 for 1-based indexing
    return compare_barcode(seq, barcodes, bc_offsets, start, end, used_offset) + 1;
}

int ShiftAdd(int sum, int digit)
{
    return sum*10 + digit;
}

//intended for umi
uint32_t get_sequence(bam1_t* b, int start, int end, int used_offset) {
    int seqlen = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
    uint8_t seq[seqlen];
    uint8_t qual[seqlen];
    bam_get_seq_qual(b, seq, qual);

    int local_start = start + used_offset;
    int local_end = end + used_offset;

    int d = 0;
    // occurrs if sequencing is truncated at beginning
    if (local_start < 0) {
        local_start = 0;
        d = abs(used_offset);
    }

    int min = 100000;
    for (int j = local_start; j <= local_end ; ++j) {
        if (int(qual[j]) < min) min = int(qual[j]);
    }

    if (min < 20) {
        return 0;
    } else {
        vector<uint8_t> umi_int(local_end - local_start + 1 + d, 0);
        size_t umi_int_size = umi_int.size();
        for (unsigned int j = 0; j < (umi_int_size - d); ++j) {
             //umi_int[j + d] = (uint8_t)bam_seqi(seq, j + local_start);
             umi_int[j + d] = seq[j + local_start];
        }
        // if offset truncates UMI, randomly assign base to truncated positions
        if (d > 0) {
            uint8_t nuc_set[] = {1,2,4,8};
            for (int i = 0; i < d; ++i) {
                umi_int[i] = nuc_set[(int)(rand() % (sizeof(nuc_set) / sizeof(int)))];
            }
        }
        uint32_t umi(accumulate(umi_int.begin(), umi_int.end(), 0, ShiftAdd));
        return umi;
    }
}

//intended for umi
string get_sequence(char* seq, char* qual, int start, int end, int min_qual, int used_offset)
{
    string sseq = string(seq);
    int local_start = start + used_offset;
    int local_end = end + used_offset;

    //int d = 0;
    // occurrs if sequencing is truncated at beginning
    if (local_start < 0)
    {
        local_start = 0;
    //    d = abs(used_offset);
    }

    int min = 100000;
    int qual_int = 0;
    for (int j = local_start; j <= local_end ; ++j)
    {
        qual_int = int(qual[j]) - 33;
        //cout << qual[j] << " " << qual_int << endl;
        if (qual_int < min_qual) return "NNNNNNNN";
    }

    // if (min < min_qual)
    // {
    //     return "NNNNNNNN";
    //} else
    //{
        return sseq.substr(local_start, local_end - local_start + 1);
    //}
}

void print_uint8 (uint8_t* arr, int seqlen, bool convert)
{
    for (int i = 0 ; i < seqlen ; ++i) {
        if (convert) {
            printf("%d", bam_seqi(arr, i));
        } else {
            printf("%d", arr[i]);
        }
    }
    printf("\n");
}

const char* convert_to_cstr(const string & s)
{
    return s.c_str();
}
