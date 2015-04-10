#include <timing.h>

void print_time(clock_t start) {
     double duration;
     duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
     cout << "[RUNTIME] " << duration << "s" << endl;
}
