/**
 * Radiation :: Output :: SoVVolume
 *
 * Welcome message and info for the 'SoVVolume' output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "SOFT.h"
#include "Tools/Radiation/Output/SoVVolume.h"

using namespace __Radiation;
using namespace std;

/**
 * Print info/welcome message.
 */
void SoVVolume::Welcome(const string& prefix) {
    slibreal_t vsize = this->arraysize;
    int unit = 0;
    char vsize_prefix[8][4] = {
        "B", "kiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB"
    };

    while (vsize > 1024.0) {
        vsize /= 1024.0;
        unit++;
    }

    SOFT::PrintInfo(prefix+"Volume array size:   %.1f %s (%zu bytes)",
        vsize, vsize_prefix[unit], this->arraysize
    );
}

