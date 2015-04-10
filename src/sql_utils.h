#ifndef SQL_UTILS_H
#define SQL_UTILS_H

#include <bam_sql_serial2.h>
#include <boost/filesystem.hpp>
#include <fstream>
#include <string>

using namespace std;

string get_selfpath();
const char* read_sql(const char* fname);

#endif
