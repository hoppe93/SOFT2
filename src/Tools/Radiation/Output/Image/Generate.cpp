/**
 * Radiation :: Output :: Image
 *
 * Generate an output image file.
 */

#include <omp.h>
#include <string>
#include <softlib/config.h>
#include <softlib/ImageGenerator/ImageGenerator.h>
#include "config.h"
#include "Tools/Radiation/Output/Image.h"

#ifdef WITH_MPI
#   include "SMPI.h"
#endif

using namespace __Radiation;
using namespace std;

slibreal_t *Image::global_image = nullptr,
           *Image::global_imageQ = nullptr,
           *Image::global_imageU = nullptr,
           *Image::global_imageV = nullptr;

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
            
            if (this->MeasuresPolarization()) {
                global_imageQ = this->imageQ;
                global_imageU = this->imageU;
                global_imageV = this->imageV;
            }
        } else {
            for (int i = 0; i < this->ntotpixels; i++)
                global_image[i] += this->image[i];

            if (this->MeasuresPolarization()) {
                for (int i = 0; i < this->ntotpixels; i++) {
                    global_imageQ[i] += this->imageQ[i];
                    global_imageU[i] += this->imageU[i];
                    global_imageV[i] += this->imageV[i];
                }
            }
            delete [] this->image;
        }
    }
}

/**
 * Generate the output file.
 * Called on the root thread only.
 */
void Image::Generate() {
#ifdef WITH_MPI
    int nprocesses, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    void *inbuf = global_image;

    if (mpi_rank == MPI_ROOT_PROCESS)
        inbuf = MPI_IN_PLACE;

    SOFT::PrintMPI("Reducing image '%s'...", this->GetName().c_str());

    unsigned int pxls = this->nrowpixels*this->ncolpixels;
    MPI_Reduce(inbuf, global_image, pxls, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD);

    if (this->MeasuresPolarization()) {
        void *inbuf_Q=global_imageQ, *inbuf_U=global_imageU, *inbuf_V=global_imageV;
        if (mpi_rank == MPI_ROOT_PROCESS) {
            inbuf_Q = MPI_IN_PLACE;
            inbuf_U = MPI_IN_PLACE;
            inbuf_V = MPI_IN_PLACE;
        }

        MPI_Reduce(inbuf_Q, global_imageQ, pxls, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD);
        MPI_Reduce(inbuf_U, global_imageU, pxls, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD);
        MPI_Reduce(inbuf_V, global_imageV, pxls, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD);
    }

    SOFT::PrintMPI("Image '%s' reduced.", this->GetName().c_str());

    if (mpi_rank != MPI_ROOT_PROCESS)
        return;
#endif

    // Write file
    string ext = this->output.substr(this->output.size()-3);

#ifdef HAS_LIBPNG
    if (ext == "png") {
        NormalizeImage();
        ImageGenerator::gen(
            ImageGenerator::IMAGE_TYPE_PNG, this->output,
            global_image, this->nrowpixels, this->ncolpixels,
            ImageGenerator::COLORMAP_GERIMAP
        );
    } else
#endif
    if (ext == "ppm") {
        NormalizeImage();
        ImageGenerator::gen(
            ImageGenerator::IMAGE_TYPE_PPM, this->output,
            global_image, this->nrowpixels, this->ncolpixels,
            ImageGenerator::COLORMAP_GERIMAP
        );
    } else {
        SFile *sf = SFile::Create(this->output, SFILE_MODE_WRITE);

        slibreal_t **img = new slibreal_t*[this->nrowpixels];
        for (int i = 0; i < this->nrowpixels; i++)
            img[i] = global_image+(i*this->ncolpixels);

        if (!this->MeasuresPolarization())
            sf->WriteArray("image", img, this->nrowpixels, this->ncolpixels);
        else {
            slibreal_t **imgQ = new slibreal_t*[this->nrowpixels];
            slibreal_t **imgU = new slibreal_t*[this->nrowpixels];
            slibreal_t **imgV = new slibreal_t*[this->nrowpixels];

            for (int i = 0; i < this->nrowpixels; i++) {
                imgQ[i] = global_imageQ+(i*this->ncolpixels);
                imgU[i] = global_imageU+(i*this->ncolpixels);
                imgV[i] = global_imageV+(i*this->ncolpixels);
            }

            sf->WriteArray("StokesI", img, this->nrowpixels, this->ncolpixels);
            sf->WriteArray("StokesQ", imgQ, this->nrowpixels, this->ncolpixels);
            sf->WriteArray("StokesU", imgU, this->nrowpixels, this->ncolpixels);
            sf->WriteArray("StokesV", imgV, this->nrowpixels, this->ncolpixels);
        }

        this->WriteCommonQuantities(sf);

        sf->Close();

        delete [] img;
    }

#ifdef COLOR_TERMINAL
    // Check if image is empty
    bool empty = true;
    unsigned int N = this->nrowpixels*this->ncolpixels;
    for (unsigned int i = 0; i < N; i++) {
        if (global_image[i] != 0) {
            empty = false;
            break;
        }
    }

    if (empty)
        SOFT::PrintInfo("Wrote \x1B[1;33mempty\x1B[0m image to '%s'.", this->output.c_str());
    else
#endif
        SOFT::PrintInfo("Wrote image to '%s'.", this->output.c_str());

    delete [] global_image;
}

void Image::NormalizeImage() {
    const unsigned int N = this->nrowpixels * this->ncolpixels;

    // Locate maximum of image
    slibreal_t imax = 0.;
    for (unsigned int i = 0; i < N; i++) {
        if (global_image[i] > imax)
            imax = global_image[i];
    }

    if (imax == 0) return;

    // Normalize all pixels to the max value
    for (unsigned int i = 0; i < N; i++)
        global_image[i] /= imax;
}

