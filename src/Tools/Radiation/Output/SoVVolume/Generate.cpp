/**
 * Radiation :: Output :: SoVVolume
 *
 * Generate an output image file.
 */

#include <omp.h>
#include <softlib/config.h>
#include "Tools/Radiation/Output/SoVVolume.h"

using namespace __Radiation;

slibreal_t *SoVVolume::global_volumearray = nullptr;

/**
 * Finish gathering output and merge output
 * stored on individual threads to the global
 * image.
 */
void SoVVolume::Finish() {
    #pragma omp critical (SoVVolume_Finish)
    {
        if (global_volumearray == nullptr) {
            global_volumearray = this->volumearray;
        } else {
            for (size_t i = 0; i < this->narrayelements; i++)
                global_volumearray[i] += this->volumearray[i];

            delete [] this->volumearray;
        }
    }
}

/**
 * Generate the output file.
 * Called on the root thread only.
 */
void SoVVolume::Generate() {
    SFile *sf = SFile::Create(this->output, SFILE_MODE_WRITE);

    sf->WriteList("volumearray", global_volumearray, this->narrayelements);
    sf->WriteList("param1", this->p1grid, this->np1);
    sf->WriteList("param2", this->p2grid, this->np2);
    sf->WriteString("param1name", Particle::GetCoordinateName(this->p1type));
    sf->WriteString("param2name", Particle::GetCoordinateName(this->p2type));

    sf->Close();

    // Check if image is empty
    bool empty = true;
    for (size_t i = 0; i < this->narrayelements; i++) {
        if (global_volumearray[i] != 0) {
            empty = false;
            break;
        }
    }

    if (empty)
#ifdef COLOR_TERMINAL
        SOFT::PrintInfo("Wrote \x1B[1;33mempty\x1B[0m volume array to '%s'.", this->output.c_str());
#else
        SOFT::PrintInfo("Wrote empty volume array to '%s'.", this->output.c_str());
#endif
    else
        SOFT::PrintInfo("Wrote volume array to '%s'.", this->output.c_str());

    delete [] global_volumearray;
}

