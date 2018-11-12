/**
 * Radiation :: Output :: Topview
 *
 * Generate a topview output file.
 */

#include <omp.h>
#include <softlib/config.h>
#include "Tools/Radiation/Output/Topview.h"

using namespace __Radiation;

slibreal_t *Topview::global_topview = nullptr;

/**
 * Finish gathering output and merge output
 * stored on individual threads to the global
 * topview.
 */
void Topview::Finish() {
    #pragma omp critical (Topview_Finish)
    {
        if (global_topview == nullptr) {
            global_topview = this->topview;
        } else {
            for (unsigned int i = 0; i < this->ntotpixels; i++)
                global_topview[i] += this->topview[i];

            delete [] this->topview;
        }
    }
}

/**
 * Generate the output file.
 * Called on the root thread only.
 */
void Topview::Generate() {
    SFile *sf = SFile::Create(this->output, SFILE_MODE_WRITE);
    slibreal_t detpos[3], detdir[3], detvisang;

    detector->GetPosition().ToArray(detpos);
    detector->GetDirection().ToArray(detdir);
    detvisang = 2.0 * detector->GetVisionAngleFOV();

    slibreal_t **img = new slibreal_t*[this->npixels];
    for (unsigned int i = 0; i < this->npixels; i++)
        img[i] = global_topview+(i*this->npixels);

    sf->WriteArray("image", img, this->npixels, this->npixels);
    sf->WriteList("detectorPosition", detpos, 3);
    sf->WriteList("detectorDirection", detdir, 3);
    sf->WriteList("detectorVisang", &detvisang, 1);

    // Include wall/separatrix data?
    slibreal_t *domain[2];
    domain[0] = magfield->GetRDomain();
    domain[1] = magfield->GetZDomain();

    sf->WriteArray("wall", domain, 2, magfield->GetNDomain());

    sf->Close();

#ifdef COLOR_TERMINAL
    // Check if image is empty
    bool empty = true;
    for (unsigned int i = 0; i < this->npixels; i++) {
        for (unsigned int j = 0; j < this->npixels; j++) {
            if (img[i][j] != 0) {
                empty = false;
                break;
            }
        }
    }

    if (empty)
        SOFT::PrintInfo("Wrote \e[1;33mempty\e[0m topview to '%s'.", this->output.c_str());
    else
#endif
    SOFT::PrintInfo("Wrote topview to '%s'.", this->output.c_str());

    delete [] img;
    delete [] global_topview;
}

