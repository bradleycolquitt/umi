#include <htslib/sam.h>
#include <bam_utils.h>
#include <iostream>
using namespace std;

int main(int argc, char** argv) {
    const char* bam_fname = argv[1];
    samFile* bam = sam_open(bam_fname, "rb");
    bam_hdr_t* header = sam_hdr_read(bam);
    hts_idx_t* idx = bam_index_load(bam_fname);

    cout << count_bam_records(idx, header) << endl;
    // for (int i=0; i < header->n_targets; ++i) {
    //         cout << header->target_len[i] << endl;
    // }


    return 0;
}
