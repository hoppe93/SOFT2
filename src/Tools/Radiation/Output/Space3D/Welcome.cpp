/**
 * Radiation :: Output :: Space3D
 *
 * Welcome message and info for the 'Space3D' output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "SOFT.h"
#include "Tools/Radiation/Output/Space3D.h"

using namespace __Radiation;
using namespace std;

/**
 * Print info/welcome message.
 */
void Space3D::Welcome(const string &prefix) {
    slibreal_t isize = this->imagesize;
    int unit = 0;
    char isize_prefix[8][4] = {
        "B", "kiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB"
    };

    while (isize > 1024.0) {
        isize /= 1024.0;
        unit++;
    }

    SOFT::PrintInfo(prefix+"S3D size:            %.1f %s (%zu bytes)",
        isize, isize_prefix[unit], this->imagesize
    );
}

