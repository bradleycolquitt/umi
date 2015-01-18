#include <fstream>
#include <bam_utils.h>
#include <htslib/sam.h>

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
