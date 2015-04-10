#include <sql_utils.h>

// /*
// ** Argument pStmt is a prepared SQL statement. This function compiles
// ** an EXPLAIN QUERY PLAN command to report on the prepared statement,
// ** and prints the report to stdout using printf().
// */
// int printExplainQueryPlan(sqlite3_stmt *pStmt){
//   const char *zSql;               /* Input SQL */
//   char *zExplain;                 /* SQL with EXPLAIN QUERY PLAN prepended */
//   sqlite3_stmt *pExplain;         /* Compiled EXPLAIN QUERY PLAN command */
//   int rc;                         /* Return code from sqlite3_prepare_v2() */

//   zSql = sqlite3_sql(pStmt);
//   if( zSql==0 ) return SQLITE_ERROR;

//   zExplain = sqlite3_mprintf("EXPLAIN QUERY PLAN %s", zSql);
//   if( zExplain==0 ) return SQLITE_NOMEM;

//   rc = sqlite3_prepare_v2(sqlite3_db_handle(pStmt), zExplain, -1, &pExplain, 0);
//   sqlite3_free(zExplain);
//   if( rc!=SQLITE_OK ) return rc;

//   while( SQLITE_ROW==sqlite3_step(pExplain) ){
//     int iSelectid = sqlite3_column_int(pExplain, 0);
//     int iOrder = sqlite3_column_int(pExplain, 1);
//     int iFrom = sqlite3_column_int(pExplain, 2);
//     const char *zDetail = (const char *)sqlite3_column_text(pExplain, 3);

//     printf("%d %d %d %s\n", iSelectid, iOrder, iFrom, zDetail);
//   }

//   return sqlite3_finalize(pExplain);
// }

string get_selfpath() {
    char buff[BUFFER_SIZE];
    ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
    if (len != -1) {
        buff[len] = '\0';
        return string(buff);
    } else {
        cerr << "File not found." << endl;
        return NULL;
    }
}

const char* read_sql(const char* fname) {
    boost::filesystem::path p(get_selfpath());
    boost::filesystem::path dir = p.parent_path();
    string path = dir.string();
    path = path + fname;

    ifstream sql_file(path);
    string sql_contents((istreambuf_iterator<char>(sql_file)), istreambuf_iterator<char>());
    const char* sql = sql_contents.c_str();
    return(sql);
    }
