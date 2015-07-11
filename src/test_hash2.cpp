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
    const char* bam_fname = "/media/data/bam/150709/ercc_cat_space/lib44-exo/lib44-exo_merge.bam";
    //const char* anno_fname = "/home/brad/src/umi/test_files/load2db_test_qsort_1M.bam.featureCounts";
    const char* fastq_fname = "/media/data/miseq/150709/lib44-exo_S1_L001_R1_001.fastq.gz";
    const char* anno_fname = "/media/data/bam/150709/ercc_cat_space/lib44-exo/lib44-exoAligned.bam.featureCounts";
    const char* dest_fname = "/home/brad/src/umi/test_files/out.txt";

    BamHash* bamhash = new BamHash(bam_fname, fastq_fname, anno_fname, dest_fname, "bc8", 8, 10, 1, 1, true);


    #ifdef PROFILE
    ProfilerStart("/home/brad/src/umi/profiling/test_hash2_1M.prof");
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
