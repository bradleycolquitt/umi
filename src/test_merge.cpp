#include <htslib/sam.h>
#include <bam_sql_serial2.h>
#include <bam_utils.h>
#include <iostream>
#include <gperftools/profiler.h>
//#include <gperftools/heap-profiler.h>


using namespace std;

int main() {
    const char* bam_fname = "/home/brad/src/umi/test_files/load2db_test_qsort_100K.bam";
    const char* db_fname = "/home/brad/src/umi/test_files/test_pe.db";
    remove(db_fname);

    BamDB* bamdb = new BamDB(bam_fname, db_fname, "bc8", 8, 17, true);
    ProfilerStart("/home/brad/src/umi/profiling/test_pe_cpu1.prof");
    //HeapProfilerStart("/home/brad/src/umi/profiling/test-heap");
    #ifdef DEBUG
    cerr << "finished align" << endl;
    #endif
    //bamdb->create_align_table();
    #ifdef DEBUG
    cerr << "finished align" << endl;
    #endif
    fill_db(bamdb);
    bamdb->create_reftable();
    ProfilerStop();
    delete bamdb;
    //HeapProfilerStop();
    return 0;
}
