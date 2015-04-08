#ifndef SQLITE_WRAPPER_H
#define SQLITE_WRAPPER_H

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
    public:
        sql_exception(int const result, char const * text)
            : code { result }
            ,  message { text }
            {}
        virtual const char* what() const throw() {
            string str = boost::str(boost::format("SQL error (%d): %s") % code % message);
            return (str.c_str());
            //return (printf("SQL error (%d): %s", code, message.c_str()));
        }
};

void execute(sqlite3* conn, const char* sql);
int prepare_statment(sqlite3* conn, const char* sql, sqlite3_stmt* stmt);
int bind(sqlite3* conn, sqlite3_stmt* stmt, int index, int value);
int bind(sqlite3* conn, sqlite3_stmt* stmt, int index, const char* text);
bool step(sqlite3* conn, sqlite3_stmt* stmt);

#endif
