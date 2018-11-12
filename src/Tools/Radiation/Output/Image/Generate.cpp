/**
 * Radiation :: Output :: Image
 *
 * Generate an output image file.
 */

#include <omp.h>
#include <softlib/config.h>
#include "Tools/Radiation/Output/Image.h"

using namespace __Radiation;

slibreal_t *Image::global_image = nullptr;

/**
 * Finish gathering output and merge output
 * stored on individual threads to the global
 * image.
 */
void Image::Finish() {
    #pragma omp critical (Image_Finish)
    {
        if (global_image == nullptr) {
            global_image = this->image;
        } else {
            for (int i = 0; i < this->ntotpixels; i++)
                global_image[i] += this->image[i];

            delete [] this->image;
        }
    }
}

/**
 * Generate the output file.
 * Called on the root thread only.
 */
void Image::Generate() {
    SFile *sf = SFile::Create(this->output, SFILE_MODE_WRITE);
    slibreal_t detpos[3], detdir[3], detvisang;

    detector->GetPosition().ToArray(detpos);
    detector->GetDirection().ToArray(detdir);
    detvisang = 2.0 * detector->GetVisionAngleFOV();

    slibreal_t **img = new slibreal_t*[this->nrowpixels];
    for (int i = 0; i < this->nrowpixels; i++)
        img[i] = global_image+(i*this->ncolpixels);

    sf->WriteArray("image", img, this->nrowpixels, this->ncolpixels);
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
    for (int i = 0; i < this->nrowpixels; i++) {
        for (int j = 0; j < this->ncolpixels; j++) {
            if (img[i][j] != 0) {
                empty = false;
                break;
            }
        }
    }

    if (empty)
        SOFT::PrintInfo("Wrote \e[1;33mempty\e[0m image to '%s'.", this->output.c_str());
    else
#endif
        SOFT::PrintInfo("Wrote image to '%s'.", this->output.c_str());

    delete [] img;
    delete [] global_image;
}

