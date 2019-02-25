/**
 * Radiation :: Output :: Green
 *
 * Welcome message and info for the 'Green' output module.
 */

#include <omp.h>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "SOFT.h"
#include "Tools/Radiation/Output/Green.h"

using namespace __Radiation;
using namespace std;

string Green::TranslateFormat(const string& fmt) {
    string s = "";

    auto capitalize = [](const string& ts) {
        string s(ts);
        for (auto& x : s)
            x = toupper(x);

        return s;
    };

    for (string::const_iterator it = fmt.begin(); it != fmt.end(); ++it) {
        if (it != fmt.begin())
            s += " x ";

        switch (*it) {
            case '1':
                s += capitalize(Particle::GetCoordinateName(this->p1type));
                break;
            case '2':
                s += capitalize(Particle::GetCoordinateName(this->p2type));
                break;
            case 'i': s += "PIXEL-I"; break;
            case 'j': s += "PIXEL-J"; break;
            case 'r': s += "RADIUS"; break;
            case 'w': s += "WAVELENGTH"; break;
            default:
                s += "<UNKN>";
                break;
        }
    }
    
    return s;
}

/**
 * Print info/welcome message.
 */
void Green::Welcome(const string &prefix) {
    size_t gf_totsize = this->fsize * omp_get_num_threads();
    slibreal_t gfsize = this->fsize;
    slibreal_t gfmem  = gf_totsize;

    int unit = 0, unitmem = 0;
    char gfsize_prefix[8][4] = {
        "B", "kiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB"
    };

    while (gfsize > 1024.0) {
        gfsize /= 1024.0;
        unit++;
    }

    while (gfmem > 1024.0) {
        gfmem /= 1024.0;
        unitmem++;
    }

    string fmt = TranslateFormat(format);
    SOFT::PrintInfo(prefix+"Function format:     %s", fmt.c_str());
    SOFT::PrintInfo(prefix+"Function size:       %.1f %s (%zu bytes)",
        gfsize, gfsize_prefix[unit], this->fsize
    );
    SOFT::PrintInfo(prefix+"Required memory:     %.1f %s (%zu bytes, %d threads)",
        gfmem, gfsize_prefix[unitmem], gf_totsize, omp_get_num_threads()
    );
}

