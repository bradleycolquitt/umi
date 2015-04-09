#include <sqlite_wrapper.h>

void execute(sqlite3* conn, const char* sql) {
    //char* err_msg = 0;
    int result = 0;
    if ((result = sqlite3_exec(conn, sql, NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception(result, sqlite3_errmsg(conn));
    }
    //sqlite3_free(err_msg);
    return;
}

// DOESN'T WORK
int prepare_statment(sqlite3* conn, const char* sql, sqlite3_stmt* stmt) {
    int result = 0;
    DEBUG_LOG(sql);
    if ((result = sqlite3_prepare_v2(conn, sql, -1, &stmt, NULL)) != SQLITE_OK) {
        throw sql_exception(result, sqlite3_errmsg(conn));
    }
    return(result);
}

int bind(sqlite3* conn, sqlite3_stmt* stmt, int index, int value) {
    int result = 0;
    if ((result = sqlite3_bind_int(stmt, index, value)) != SQLITE_OK) {
        throw sql_exception(result, sqlite3_errmsg(conn));
    }
    return result;
}

// DOESN'T WORK
int bind(sqlite3* conn, sqlite3_stmt* stmt, int index, const char* text) {
    int result = 0;
    cerr << text << endl;
    DEBUG_LOG(text);
    if ((result = sqlite3_bind_text(stmt, index, text, -1, SQLITE_TRANSIENT)) != SQLITE_OK) {
        throw sql_exception(result, sqlite3_errmsg(conn));
    }
    return result;
}

bool step(sqlite3* conn, sqlite3_stmt* stmt) {
    auto const result = sqlite3_step(stmt);
    if (result == SQLITE_ROW) return true;
    if (result == SQLITE_DONE) return false;

    throw sql_exception { result, sqlite3_errmsg(conn)};
}
