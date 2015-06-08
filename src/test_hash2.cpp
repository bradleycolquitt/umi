#include <htslib/sam.h>
#include <bam_hash2.h>
#include <bamrecord.h>
#include <bam_utils.h>
#include <timing.h>
#include <iostream>
#include <fstream>
#include <ctime>

//#define PROFILE
#ifdef PROFILE
#include <gperftools/profiler.h>
#include <gperftools/heap-profiler.h>
#endif


using namespace std;

int main() {
    //const char* bam_fname = "/home/brad/src/umi/test_files/load2db_test_qsort_1M.bam";
    const char* bam_fname = "/media/data/bam/150403/ercc_cat_space/lib40/lib40_clip18Aligned.bam";
    //const char* anno_fname = "/home/brad/src/umi/test_files/load2db_test_qsort_1M.bam.featureCounts";
    const char* anno_fname = "/media/data/bam/150403/ercc_cat_space/lib40/lib40_clip18Aligned.bam.featureCounts";
    const char* dest_fname = "/home/brad/src/umi/test_files/out.txt";

    BamHash* bamhash = new BamHash(bam_fname, anno_fname, dest_fname, "bc8", 8, 10, 1, 1);


    #ifdef PROFILE
    ProfilerStart("/home/brad/src/umi/profiling/test_hash_1M.prof");
    //HeapProfilerStart("/home/brad/src/umi/profiling/test_merge.heap");
    #endif

    clock_t start = clock();
    hash_annotation(bamhash);
    print_time(start);

    bamhash->hash_reads();
    print_time(start);

    bamhash->print_results();
    print_time(start);

    #ifdef PROFILE
    ProfilerStop();
    //HeapProfilerStop();
    #endif
    delete bamhash;

    return 0;
}
