#include <sqlite_wrapper.h>
#include <thread_pool.h>
#include <sqlite3.h>
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/exception_ptr.hpp>
#include <fstream>
#include <mutex>

class FStreamWriter
{
    private:
        std::fstream *f;
        std::mutex mtx;
    public:
        FStreamWriter(std::fstream *f) : f(f) {}
        template <typename T>
        FStreamWriter &operator<<(const T &x)
        {
            mtx.lock();
            //std::lock_guard<std::mutex> lock(mtx); // <- Mutex to make it safe-thread
            (*f) << x;
            mtx.unlock();
            return *this;
        }
        void flush() {
            f->flush();
        }
};

/*!
Execute SQL function
 @param tid
 @param sql
 @param db_fname
*/
void execute_sql(int tid, const char* sql, const char* db_fname);

/*

*/
