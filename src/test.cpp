#include <htslib/sam.h>
#include <bam_sql.h>
#include <bam_utils.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    const char* bam_fname = argv[1];
    const char* dest_fname = "test.db";
    samFile* bam = sam_open(bam_fname, "rb");
    bam_hdr_t* header = sam_hdr_read(bam);
    hts_idx_t* idx = bam_index_load(bam_fname);
    cout << count_bam_records(idx, header) << endl;
    BamDB* bamdb = new BamDB(bam_fname, dest_fname, "bc8");
    create_table(bamdb);
    fill_db(bamdb);
    return 0;
}
