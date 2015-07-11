#include <sqlite_wrapper.h>

void open_connection(const char* fname, sqlite3** conn, bool for_speed) {
    if (sqlite3_open(fname, conn) != SQLITE_OK) {
        throw sql_exception(-1, sqlite3_errmsg(*conn));
    }

    if (for_speed) {
        execute(*conn, "PRAGMA journal_mode = MEMORY");
        execute(*conn, "PRAGMA synchronous = OFF");
    }
}

void execute(sqlite3* conn, const char* sql) {
    //char* err_msg = 0;
    int result = 0;
    if ((result = sqlite3_exec(conn, sql, NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception(result, sqlite3_errmsg(conn));
    }
    //sqlite3_free(err_msg);
    return;
}

void execute(sqlite3** conn, const char* sql) {
    //char* err_msg = 0;
    int result = 0;
    if ((result = sqlite3_exec(*conn, sql, NULL, NULL, NULL)) != SQLITE_OK) {
        throw sql_exception(result, sqlite3_errmsg(*conn));
    }
    //sqlite3_free(err_msg);
    return;
}

// DOESN'T WORK
int prepare_statment(sqlite3* conn, const char* sql, sqlite3_stmt** stmt) {
    int result = 0;
    //DEBUG_LOG(sql);
    if ((result = sqlite3_prepare_v2(conn, sql, -1, stmt, NULL)) != SQLITE_OK) {
        throw sql_exception(result, sqlite3_errmsg(conn));
    }
    return(result);
}

int prepare_statment(sqlite3** conn, const char* sql, sqlite3_stmt** stmt, const char** tail) {
    int result = 0;
    //DEBUG_LOG(sql);
    if ((result = sqlite3_prepare_v2(*conn, sql, -1, stmt, tail)) != SQLITE_OK) {
        throw sql_exception(result, sqlite3_errmsg(*conn), "prepare_statement");
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
    //DEBUG_LOG(text);
    if ((result = sqlite3_bind_text(stmt, index, text, -1, SQLITE_TRANSIENT)) != SQLITE_OK) {
        throw sql_exception(result, sqlite3_errmsg(conn));
    }
    return result;
}

bool step(sqlite3 *& conn, sqlite3_stmt** stmt) {
    auto const result = sqlite3_step(*stmt);
    if (result == SQLITE_ROW) return true;
    if (result == SQLITE_DONE) return false;

    throw sql_exception { result, sqlite3_errmsg(conn), "step"};
}

bool step_multiple(sqlite3* conn, const char* sql) {
    sqlite3_stmt* stmt;
    int result = 0;
    const char* tail;

    if ((result = sqlite3_prepare_v2(conn, sql, -1, &stmt, &tail)) != SQLITE_OK) {
        throw sql_exception { result, sqlite3_errmsg(conn)};
    }

    if ((result = sqlite3_step(stmt)) <100 ) {
        throw sql_exception { result, sqlite3_errmsg(conn)};
    }

    sqlite3_reset(stmt);
    while (*tail != '\0') {
        //DEBUG_LOG(tail);
        result = sqlite3_prepare_v2(conn, tail, -1, &stmt, &tail);
        if (result != SQLITE_OK) {
            throw sql_exception { result, sqlite3_errmsg(conn)};
        }
        sqlite3_step(stmt);
        sqlite3_reset(stmt);
    }
    sqlite3_finalize(stmt);
}

void reset(sqlite3** conn, sqlite3_stmt** stmt_p) {
    int result = 0;
        if ((result = sqlite3_reset(*stmt_p)) != SQLITE_OK) {
        throw sql_exception(result, sqlite3_errmsg(*conn), "reset");
    }
}
