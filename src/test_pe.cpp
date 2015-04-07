#include <htslib/sam.h>
#include <bam_sql.h>
#include <bam_utils.h>
#include <iostream>
#include <gperftools/profiler.h>
//#include <gperftools/heap-profiler.h>


using namespace std;

int main() {
    const char* bam_fname = "/home/brad/src/umi/test_files/load2db_test_qsort.bam";
    const char* db_fname = "/home/brad/src/umi/test_files/test_pe.db";
    remove(db_fname);

    // samFile* bam = sam_open(bam_fname, "rb");
    // bam_hdr_t* header = sam_hdr_read(bam);
    // hts_idx_t* idx = bam_index_load(bam_fname);
    // cout << count_bam_records(idx, header) << endl;
    BamDB* bamdb = new BamDB(bam_fname, db_fname, "bc8", 8, 17, true);
    //ProfilerStart("/home/brad/src/umi/profiling/test_dna3_cpu.prof");
    //HeapProfilerStart("/home/brad/src/umi/profiling/test-heap");
    #ifdef DEBUG
    cerr << "finished align" << endl;
    #endif
    bamdb->create_align_table();
    #ifdef DEBUG
    cerr << "finished align" << endl;
    #endif
    fill_db(bamdb);
    bamdb->create_reftable();
    //ProfilerStop();
    delete bamdb;
    //HeapProfilerStop();
    return 0;
}
