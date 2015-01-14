#include <htslib/sam.h>

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
    vector<string>::iterator it = s.begin();
    vector<string>::iterator it_end = s.end();
    vector<int> out(s.size());
    for (int i=0; i < s.size() ; ++i) {
        switch ( s[i])
        {
            case 'A': out[i] = 1;
            case 'C': out[i] = 2;
            case 'G': out[i] = 4;
            case 'T': out[i] = 8;
            case 'N': out[i] = 15;
        }
    }
    return out;
}
