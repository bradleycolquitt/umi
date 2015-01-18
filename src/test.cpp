#include <htslib/sam.h>
#include <bam_sql.h>
#include <bam_utils.h>
#include <iostream>
#include <gperftools/profiler.h>
//#include <gperftools/heap-profiler.h>


using namespace std;

int main(int argc, char** argv) {
    const char* bam_fname = argv[1];
    const char* dest_fname = "test.db";

    remove(dest_fname);

    samFile* bam = sam_open(bam_fname, "rb");
    bam_hdr_t* header = sam_hdr_read(bam);
    hts_idx_t* idx = bam_index_load(bam_fname);
    cout << count_bam_records(idx, header) << endl;
    BamDB* bamdb = new BamDB(bam_fname, dest_fname, "bc8");
    ProfilerStart("/home/brad/src/umi/profiling/test_dna3_cpu.prof");
    //HeapProfilerStart("/home/brad/src/umi/profiling/test-heap");
    create_table(bamdb);
    fill_db(bamdb);
    bamdb->create_reftable();
    ProfilerStop();
    delete bamdb;
    //HeapProfilerStop();
    return 0;
}
