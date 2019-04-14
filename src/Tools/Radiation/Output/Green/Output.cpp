/**
 * Radiation :: Output :: Green
 *
 * Generate an output Green's function.
 */

#include <omp.h>
#include <softlib/config.h>
#include "config.h"
#include "Tools/Radiation/Output/Green.h"

#ifdef WITH_MPI
#   include <mpi.h>
#   include "SMPI.h"
#endif

using namespace __Radiation;
using namespace std;

slibreal_t *Green::global_function = nullptr;

/**
 * Finish gathering output and merge output
 * stored on individual threads to the global
 * Green's function.
 */
void Green::Finish() {
    #pragma omp critical (Green_Finish)
    {
        if (global_function == nullptr) {
            global_function = this->function;
        } else {
            for (unsigned int i = 0; i < this->fsize; i++)
                global_function[i] += this->function[i];

            delete [] this->function;
        }
    }
}

/**
 * Generate the output file.
 * Called on the root thread only.
 */
void Green::Generate() {
#ifdef WITH_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    string output_name = this->output;

    if (this->mpi_output_mode == MPI_Output_Mode::CHUNKED) {
        output_name = GetChunkedName(mpi_rank);

        SOFT::PrintMPI("Saving Green's function chunk to '%s'...", output_name.c_str());
    } else {
        SOFT::PrintMPI("Reducing Green's function '%s'...", this->GetName().c_str());

        void *inbuf = global_function;
        if (mpi_rank == MPI_ROOT_PROCESS)
            inbuf = MPI_IN_PLACE;

        MPI_Reduce(
            inbuf, global_function, this->fsize, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD
        );

        SOFT::PrintMPI("Green's function '%s' reduced.", this->GetName().c_str());

        if (mpi_rank != MPI_ROOT_PROCESS)
            return;
    }
#else
    string output_name = this->output;
#endif
    SFile *sf = SFile::Create(output_name, SFILE_MODE_WRITE);

    // Store as linear or multi-dimensional array?
    if (this->storeFAsLinearArray)
        sf->WriteList("func", global_function, fsize);
    else
        sf->WriteMultiArray("func", global_function, this->ndimensions, this->dimensions);

    sf->WriteString("type", this->format);
    sf->WriteList("wavelengths", this->detector->GetWavelengths(), this->nw);

    if (this->hasI)
        sf->WriteScalar("rowpixels", this->nrowpixels);
    if (this->hasJ)
        sf->WriteScalar("colpixels", this->ncolpixels);

    if (this->storeStokesParameters)
        sf->WriteScalar("stokesparams", 1.0);
	else
		sf->WriteScalar("stokesparams", 0.0);

	this->WriteCommonQuantities(sf);

    sf->Close();

    SOFT::PrintInfo("Wrote Green's function to '%s'.", this->output.c_str());

    delete [] global_function;
}

/**
 * Returns the actual file name to use when in
 * MPI_Output_Mode "CHUNKED". In this case, each MPI process
 * generates its own output file, with the process rank
 * inserted in the file name.
 * 
 * rank: MPI process rank.
 */
string Green::GetChunkedName(const int rank) {
    string::size_type n = this->output.find('.');
    string newName(this->output), srank = to_string(rank);

    newName.insert(n, srank);

    return newName;
}

