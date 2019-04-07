/**
 * Radiation :: Output :: SoVVolume
 *
 * Generate an output image file.
 */

#include <cstdlib>
#include <omp.h>
#include <softlib/config.h>
#include "Tools/Radiation/Output/SoVVolume.h"

#ifdef WITH_MPI
#   include <mpi.h>
#   include "SMPI.h"
#endif

using namespace __Radiation;
using namespace std;

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
#ifdef WITH_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    void *inbuf = global_volumearray;
    if (mpi_rank == MPI_ROOT_PROCESS)
        inbuf = MPI_IN_PLACE;
    
    SOFT::PrintMPI("Reducing SoV volume '%s'...", this->GetName().c_str());
    MPI_Reduce(inbuf, global_volumearray, this->narrayelements, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD);
    SOFT::PrintMPI("SoV volume '%s' reduced.", this->GetName().c_str());

    if (mpi_rank != MPI_ROOT_PROCESS)
        return;
#endif

    SFile *sf = SFile::Create(this->output, SFILE_MODE_WRITE);

    //sf->WriteList("volumearray", global_volumearray, this->narrayelements);
    sfilesize_t ndims = 2, dims[2] = {this->np1, this->np2};
    sf->WriteMultiArray("volumearray", global_volumearray, ndims, dims);

    sf->WriteList("param1", this->p1grid, this->np1);
    sf->WriteList("param2", this->p2grid, this->np2);
    sf->WriteString("param1name", Particle::GetCoordinateName(this->p1type));
    sf->WriteString("param2name", Particle::GetCoordinateName(this->p2type));

	this->WriteCommonQuantities(sf);

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

