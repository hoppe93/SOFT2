/**
 * Radiation :: Output :: Spectrum
 *
 * Generate an output spectrum file.
 */

#include <omp.h>
#include <softlib/config.h>
#include "Tools/Radiation/Output/Spectrum.h"

#ifdef WITH_MPI
#   include <mpi.h>
#   include "SMPI.h"
#endif

using namespace __Radiation;

slibreal_t
    *Spectrum::global_I = nullptr,
    *Spectrum::global_Q = nullptr,
    *Spectrum::global_U = nullptr,
    *Spectrum::global_V = nullptr;

/**
 * Finish gathering output and merge output
 * stored on individual threads to the global
 * spectrum.
 */
void Spectrum::Finish() {
    #pragma omp critical (Spectrum_Finish)
    {
        if (global_I == nullptr) {
            global_I = this->I;
            if (this->MeasuresPolarization()) {
                global_Q = this->Q;
                global_U = this->U;
                global_V = this->V;
            }
        } else {
            for (unsigned int i = 0; i < this->nwavelengths; i++)
                global_I[i] += this->I[i];

            delete [] this->I;

            if (this->MeasuresPolarization()) {
                for (unsigned int i = 0; i < this->nwavelengths; i++) {
                    global_Q[i] += this->Q[i];
                    global_U[i] += this->U[i];
                    global_V[i] += this->V[i];
                }

                delete [] this->Q;
                delete [] this->U;
                delete [] this->V;
            }
        }
    }
}

/**
 * Generate the output file.
 * Called on the root thread only.
 */
void Spectrum::Generate() {
#ifdef WITH_MPI
    int nprocesses, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    void *inbuf = global_I;

    if (mpi_rank == MPI_ROOT_PROCESS)
        inbuf = MPI_IN_PLACE;

    SOFT::PrintMPI("Reducing spectrum '%s'...", this->GetName().c_str());

    MPI_Reduce(inbuf, global_I, this->nwavelengths, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD);

    if (this->MeasuresPolarization()) {
        void *inbuf_Q=global_Q, *inbuf_U=global_U, *inbuf_V=global_V;
        if (mpi_rank == MPI_ROOT_PROCESS) {
            inbuf_Q = MPI_IN_PLACE;
            inbuf_U = MPI_IN_PLACE;
            inbuf_V = MPI_IN_PLACE;
        }

        MPI_Reduce(inbuf_Q, global_Q, this->nwavelengths, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD);
        MPI_Reduce(inbuf_U, global_U, this->nwavelengths, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD);
        MPI_Reduce(inbuf_V, global_V, this->nwavelengths, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD);
    }

    SOFT::PrintMPI("Spectrum '%s' reduced.", this->GetName().c_str());

    if (mpi_rank != MPI_ROOT_PROCESS)
        return;
#endif
    SFile *sf = SFile::Create(this->output, SFILE_MODE_WRITE);

    sf->WriteList("wavelengths", this->wavelengths, this->nwavelengths);
    sf->WriteList("I", global_I, this->nwavelengths);
    if (this->MeasuresPolarization()) {
        sf->WriteList("Q", global_Q, this->nwavelengths);
        sf->WriteList("U", global_U, this->nwavelengths);
        sf->WriteList("V", global_V, this->nwavelengths);
    }

	this->WriteCommonQuantities(sf);

    sf->Close();

    SOFT::PrintInfo("Wrote spectrum to '%s'.", this->output.c_str());

    if (this->MeasuresPolarization()) {
        delete [] global_V;
        delete [] global_U;
        delete [] global_Q;
    }

    delete [] global_I;
}

