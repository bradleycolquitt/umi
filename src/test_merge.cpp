#include <htslib/sam.h>
#include <dbrecord.h>
#include <bam_sql_serial2.h>
#include <bam_utils.h>
#include <iostream>
#include <sqlite3.h>
//#define PROFILE

#ifdef PROFILE
#include <gperftools/profiler.h>
#include <gperftools/heap-profiler.h>
#endif


using namespace std;

int main() {
    const char* bam_fname = "/home/brad/src/umi/test_files/load2db_test_qsort_100K.bam";
    const char* db_fname = "/home/brad/src/umi/test_files/test_merge.db";
    remove(db_fname);

    BamDB* bamdb = new BamDB(bam_fname, db_fname, "bc8", 8, 17, true);

    #ifdef PROFILE
    //ProfilerStart("/home/brad/src/umi/profiling/test_merge_1M.prof");
    HeapProfilerStart("/home/brad/src/umi/profiling/test_merge.heap");
    #endif
    #ifdef DEBUG
    cerr << "finished align" << endl;
    #endif

    #ifdef DEBUG
    cerr << "finished align" << endl;
    #endif

    // dbRecord1* record1 = new dbRecord1(bamdb);
    // record1->insert_to_db();
    //sqlite3_exec(bamdb->get_conn(), "CREATE TABLE tests (test text);", NULL, NULL, NULL);
    fill_db(bamdb);
    bamdb->create_reftable();
    #ifdef PROFILE
    //ProfilerStop();
    HeapProfilerStop();
    #endif
    delete bamdb;

    return 0;
}
