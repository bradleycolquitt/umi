#include <sqlite_wrapper.h>

// // struct statement_handle_traits
// {
//     using pointer = sqlite3_stmt *;
//     static auto invalid() noexcept {
//         return nullptr;
//     }
//     static auto close(pointer value) noexcept {
//         VERIFY_(SQLITE_OK, sqlite3_finalize(value));
//     }
// };

// using statement_handle = unique_handle<statement_handle_traits>;

void execute(sqlite3* conn, const char* sql) {
    char* err_msg;
    int result = 0;
    if ((result = sqlite3_exec(conn, sql, NULL, NULL, &err_msg)) != SQLITE_OK) {
        throw sql_exception(result, err_msg);
    }
    sqlite3_free(err_msg);
    return;
}

int prepare_statment(sqlite3* conn, const char* sql, sqlite3_stmt* stmt) {
    int result = 0;
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

int bind(sqlite3* conn, sqlite3_stmt* stmt, int index, const char* text) {
    int result = 0;
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
