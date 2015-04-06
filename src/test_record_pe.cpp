#include <htslib/sam.h>
#include <bam_sql.h>
#include <bam_utils.h>
#include <iostream>
#include <gperftools/profiler.h>
//#include <gperftools/heap-profiler.h>

using namespace std;

int main(int argc, char** argv) {
    const char* bam_fname = argv[1];
    int umi_start = atoi(argv[2]);
    int umi_end = atoi(argv[3]);
    const char* dest_fname = "test.db";

    remove(dest_fname);

    samFile* bam = sam_open(bam_fname, "rb");
    bam_hdr_t* header = sam_hdr_read(bam);
    hts_idx_t* idx = bam_index_load(bam_fname);

    BamDB* bamdb = new BamDB(bam_fname, dest_fname, "bc8", 8, 17, true);

    bamdb->get_record()->pos_head2 = 1;
    bamdb->get_record()->pos_tail2 = 50;
    bamdb->get_record()->insert = 100;

    cout << bamdb->get_record()->pos_head2 << " " << bamdb->get_record()->pos_tail2
         << " " << bamdb->get_record()->insert << endl;
    //ProfilerStart("/home/brad/src/umi/profiling/test_dna3_cpu.prof");
    //HeapProfilerStart("/home/brad/src/umi/profiling/test-heap");
    // create_table(bamdb);
    // fill_db(bamdb);
    // bamdb->create_reftable();
    // ProfilerStop();
    // delete bamdb;
    //HeapProfilerStop();
    return 0;
}
