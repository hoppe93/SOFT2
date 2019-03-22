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

    slibreal_t **img = new slibreal_t*[this->npixels];
    for (unsigned int i = 0; i < this->npixels; i++)
        img[i] = global_topview+(i*this->npixels);

    sf->WriteArray("image", img, this->npixels, this->npixels);

	this->WriteCommonQuantities(sf);

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
        SOFT::PrintInfo("Wrote \x1B[1;33mempty\x1B[0m topview to '%s'.", this->output.c_str());
    else
#endif
    SOFT::PrintInfo("Wrote topview to '%s'.", this->output.c_str());

    delete [] img;
    delete [] global_topview;
}

