#include <htslib/sam.h>
#include <dbrecord.h>
#include <bam_sql_merge.h>
#include <bam_utils.h>
#include <iostream>
#include <sqlite3.h>

#define PROFILE
#ifdef PROFILE
#include <gperftools/profiler.h>
#include <gperftools/heap-profiler.h>
#endif


using namespace std;

int main() {
    const char* bam_fname = "/home/brad/src/umi/test_files/load2db_test_qsort_1M.bam";
    const char* db_fname = "/home/brad/src/umi/test_files/test_merge.db";
    remove(db_fname);

    BamDB* bamdb = new BamDB(bam_fname, db_fname, "bc8", 8, 10, false, true);

    #ifdef PROFILE
    ProfilerStart("/home/brad/src/umi/profiling/test_merge_1M.prof");
    //HeapProfilerStart("/home/brad/src/umi/profiling/test_merge.heap");
    #endif

    fill_db(bamdb);
//    bamdb->create_reftable();
    #ifdef PROFILE
    ProfilerStop();
    //HeapProfilerStop();
    #endif
    delete bamdb;

    return 0;
}
