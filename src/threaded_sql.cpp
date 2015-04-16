#include <threaded_sql.h>

//static const char * sql = "SELECT count(*) FROM collapsed WHERE tid=? ;";

static const char * sql = "SELECT bc, gene_id, sum(unique_umi) as unique_umi, sum(total_umi) as total_umi FROM\
(\
SELECT DISTINCT collapsed.bc, anno.genes_ref.gene_id, grouped.unique_umi, grouped.total_umi FROM\
    anno.genes_ref JOIN\
      anno.gene_rtree ON anno.genes_ref.rowid = anno.gene_rtree.id\
      LEFT JOIN\
      collapsed ON\
          anno.genes_ref.tid=collapsed.tid\
          AND anno.genes_ref.strand=collapsed.strand\
      JOIN pos_rtree\
          ON collapsed.rowid = pos_rtree.id\
          AND\
          (\
               pos_rtree.lpos1 BETWEEN anno.gene_rtree.start AND anno.gene_rtree.end\
               OR pos_rtree.rpos2 BETWEEN anno.gene_rtree.start AND anno.gene_rtree.end\
          )\
       JOIN grouped ON collapsed.pos_id=grouped.pos_id\
       WHERE\
           anno.genes_ref.build='maker'\
           AND anno.genes_ref.tid=? \
) AS tmp GROUP BY bc, gene_id";

void execute_sql(int tid, const char * sql, const char * db_fname,
                 FStreamWriter& os, boost::exception_ptr& error)
{
    cout << "tid: " << tid << endl;
    sqlite3 * conn;
    try {
        open_connection(db_fname, &conn, false);
        error = boost::exception_ptr();
    } catch ( ... ) {
        error = boost::current_exception();
        cout << tid << " err opening" << endl;
    }
    cout << "here1" << endl;
    execute(conn, "ATTACH DATABASE '/media/data/db/anno.db' as anno;");

    sqlite3_stmt * stmt;
    int rc = 0;
    if ((rc = sqlite3_prepare_v2(conn, sql, -1, &stmt, NULL)) != SQLITE_OK) {
        cout << tid << " err prep " << sqlite3_errmsg(conn) << endl;
        throw sql_exception { rc, sqlite3_errmsg(conn)};
    }
    cout << "here2" << endl;

    if ((rc = sqlite3_bind_int(stmt, 1, tid)) != SQLITE_OK) {
        cout << tid << " bind err " << sqlite3_errmsg(conn) << endl;
    }

    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW || rc == SQLITE_DONE) {
        //cout << "here3" << endl;
        if (rc < 100) {
            cout << tid << " err step: " <<  sqlite3_errmsg(conn) << endl;
            throw sql_exception { rc, sqlite3_errmsg(conn)};
        } else {
            //cout << rc << endl;
            //int num_cols = sqlite3_column_count(stmt);
            //cout << num_cols << endl;
            int result = sqlite3_column_int(stmt, 1);
            os << sqlite3_column_int(stmt, 0) << "\t";
            os << sqlite3_column_text(stmt, 1) << "\t";
            os << sqlite3_column_int(stmt, 2) << "\t";
            os << sqlite3_column_int(stmt, 3) << "\t";
            os << "\n";
            os.flush();
        }
    }
    sqlite3_finalize(stmt);
    sqlite3_close(conn);
 }

int main(int argc, char ** argv) {

    /* Annotation db */
    const char * db_fname = argv[1];
    sqlite3 * conn;

    open_connection(db_fname, &conn, false);

    /* Extract number of tids from annotation db*/
    sqlite3_stmt * stmt;
    sqlite3_prepare_v2(conn, "SELECT count(*) FROM reference;", -1, &stmt, NULL);
    sqlite3_step(stmt);
    int num_tid = sqlite3_column_int(stmt, 0);
    cerr << num_tid << endl;
    sqlite3_finalize(stmt);
    sqlite3_close(conn);

    /* Set up thread pool */
    int n_threads = atoi(argv[2]);
    thread_pool pool (n_threads);
    boost::exception_ptr error;

    /* Set ofstream */
    fstream ofs("test.txt", ios::out);
    FStreamWriter fwriter(&ofs);

    /* For each each tid start thread */
    int rc = 0;
    for (int i = 0; i < 10 ; ++i) {
        //cerr << "here" << endl;
        try
        {
            if ((rc = pool.run_task(boost::bind(execute_sql, i, sql,
                                    db_fname, std::ref(fwriter), std::ref(error)))) < 0)
            {
                --i;
            }
            //cerr << "here2" << endl;
        } catch (sql_exception &e) {
            cerr << e.what() << endl;
        }
    }
    ofs.flush();
    ofs.close();
}
