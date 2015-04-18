#include <htslib/sam.h>
#include <bamrecord.h>
#include <bam_hash.h>
#include <bam_utils.h>
#include <iostream>

#define PROFILE
#ifdef PROFILE
#include <gperftools/profiler.h>
#include <gperftools/heap-profiler.h>
#endif


using namespace std;

int main() {
    const char* bam_fname = "/home/brad/src/umi/test_files/load2db_test_qsort_1M.bam";
    const char* anno_fname = "/home/brad/src/umi/test_files/load2db_test_qsort_1M.bam.featureCounts";
    const char* dest_fname = "/home/brad/src/umi/test_files/out.txt";

    // ofstream ofs("test.txt");
    // ofs << "test" << endl;
    // ofs.close();

    BamHash* bamhash = new BamHash(bam_fname, anno_fname, dest_fname, "bc8", 8, 10);


    #ifdef PROFILE
    ProfilerStart("/home/brad/src/umi/profiling/test_hash_1M.prof");
    //HeapProfilerStart("/home/brad/src/umi/profiling/test_merge.heap");
    #endif
    hash_annotation(bamhash);
    bamhash->hash_reads();
    bamhash->print_results();
    //print_results(bamhash);
    //fill_db(bamdb);
//    bamdb->create_reftable();
    #ifdef PROFILE
    ProfilerStop();
    //HeapProfilerStop();
    #endif
    delete bamhash;

    return 0;
}
