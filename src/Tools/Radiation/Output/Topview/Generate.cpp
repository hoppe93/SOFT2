/**
 * Radiation :: Output :: Topview
 *
 * Generate a topview output file.
 */

#include <omp.h>
#include <softlib/config.h>
#include "Tools/Radiation/Output/Topview.h"

#ifdef WITH_MPI
#   include "SMPI.h"
#endif

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
#ifdef WITH_MPI
    int nprocesses, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    void *inbuf = global_topview;

    if (mpi_rank == MPI_ROOT_PROCESS)
        inbuf = MPI_IN_PLACE;

    SOFT::PrintMPI("Reducing topview '%s'...", this->GetName().c_str());

    MPI_Reduce(inbuf, global_topview, this->npixels*this->npixels, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD);

    SOFT::PrintMPI("Topview '%s' reduced.", this->GetName().c_str());

    if (mpi_rank != MPI_ROOT_PROCESS)
        return;
#endif

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

