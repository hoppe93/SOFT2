/**
 * Radiation :: Output :: Green
 *
 * Welcome message and info for the 'Green' output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "SOFT.h"
#include "Tools/Radiation/Output/Green.h"

using namespace __Radiation;
using namespace std;

/**
 * Print info/welcome message.
 */
void Green::Welcome(const string &prefix) {
    slibreal_t gfsize = this->fsize;
    int unit = 0;
    char gfsize_prefix[8][4] = {
        "B", "kiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB"
    };

    while (gfsize > 1024.0) {
        gfsize /= 1024.0;
        unit++;
    }

    SOFT::PrintInfo(prefix+"Function format:     %s", format.c_str());
    SOFT::PrintInfo(prefix+"Function size:       %.1f %s (%zu bytes)",
        gfsize, gfsize_prefix[unit], this->fsize
    );
}

