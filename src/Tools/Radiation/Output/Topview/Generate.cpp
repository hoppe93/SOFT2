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
    unsigned int nwall = magfield->GetNDomain();
    slibreal_t **domain = new slibreal_t*[2], *rwall, *zwall;
    domain[0] = new slibreal_t[2*nwall];
    domain[1] = domain[0]+nwall;

    rwall = magfield->GetRDomain();
    zwall = magfield->GetZDomain();
    for (unsigned int i = 0; i < nwall; i++) {
        domain[0][i] = rwall[i];
        domain[1][i] = zwall[i];
    }

    sf->WriteArray("wall", domain, 2, nwall);

    sf->Close();

    delete [] domain[0];
    delete [] domain;

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
        SOFT::PrintInfo("Wrote \x1B[1;33mempty\x1B[0m topview to '%s'.", this->output.c_str());
    else
#endif
    SOFT::PrintInfo("Wrote topview to '%s'.", this->output.c_str());

    delete [] img;
    delete [] global_topview;
}

