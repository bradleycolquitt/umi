#ifndef SQLITE_WRAPPER_H
#define SQLITE_WRAPPER_H

#include <bam_sql_merge.h>
#include <sqlite3.h>
#include <boost/format.hpp>
#include <exception>
#include <string>

using namespace std;
//using namespace string;
//using namespace boost::format
//struct statement_handle_traits;

class sql_exception: public exception {
    private:
        int code;
        string message;
        string message2;
    public:
        sql_exception(int const result, char const * text)
            : code { result }
            ,  message { text }
            , message2 { "" }
            {}
        sql_exception(int const result, char const * text, char const * text2)
            : code { result }
            , message { text }
            , message2 { text2 }
            {}
        virtual const char* what() const throw() {
            string str = boost::str(boost::format("SQL error [%s] (%d): %s")
                         % message2 % code % message);
            return (str.c_str());
            //return (printf("SQL error (%d): %s", code, message.c_str()));
        }
};

void open_connection(const char* fname, sqlite3** conn, bool for_speed);
void execute(sqlite3* conn, const char* sql);
void execute(sqlite3** conn, const char* sql);
int prepare_statment(sqlite3* conn, const char* sql, sqlite3_stmt** stmt);
int prepare_statment(sqlite3** conn, const char* sql, sqlite3_stmt** stmt, const char** tail);
int bind(sqlite3* conn, sqlite3_stmt* stmt, int index, int value);
int bind(sqlite3* conn, sqlite3_stmt* stmt, int index, const char* text);
bool step(sqlite3 *& conn, sqlite3_stmt** stmt);
bool step_multiple(sqlite3* conn, const char* sql);
void reset(sqlite3** conn, sqlite3_stmt** stmt_p);
#endif
