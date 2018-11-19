/**
 * Print a welcome message for the 'Orbits' tool.
 */

#include <string>
#include "SOFT.h"
#include "Tools/Orbits.h"

using namespace std;

void Orbits::Welcome(const string &prefix) {
    size_t nbytes = this->allocatedBytes;
    slibreal_t ebytes = (slibreal_t)nbytes;
    int n = 0;
    const int NSUFFIX = 9;
    char suffix[NSUFFIX][4] = { "B", "kiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB" };

    while (ebytes > 1024 && n < NSUFFIX) {
        ebytes /= 1024;
        n++;
    }

    SOFT::PrintInfo(prefix+"Allocated memory:   %.1f %s (%zu bytes)", ebytes, suffix[n], nbytes);
}
