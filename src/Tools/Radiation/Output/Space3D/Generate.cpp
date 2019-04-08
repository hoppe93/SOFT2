/**
 * Radiation :: Output :: Space3D
 *
 * Generate an output S3D file.
 */

#include <omp.h>
#include <softlib/config.h>
#include "Tools/Radiation/Output/Space3D.h"

#ifdef WITH_MPI
#   include <mpi.h>
#   include "SMPI.h"
#endif

using namespace __Radiation;
using namespace std;

/**
 * Finish processing on each thread.
 * For Space3D, we don't need to do anything
 * since all threads work on the same image.
 */
void Space3D::Finish() { }

/**
 * Generate the output file.
 * Called on the root thread only.
 */
void Space3D::Generate() {
#ifdef WITH_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    void *inbuf = this->s3dimage;

    if (mpi_rank == MPI_ROOT_PROCESS)
        inbuf = MPI_IN_PLACE;

    SOFT::PrintMPI("Reducing S3D '%s'...", this->GetName().c_str());
    MPI_Reduce(inbuf, this->s3dimage, this->imagesize, SMPI::MPI_SLIBREAL_T, SMPI::SUM, MPI_ROOT_PROCESS, MPI_COMM_WORLD);
    SOFT::PrintMPI("S3D '%s' reduced.", this->GetName().c_str());

    if (mpi_rank != MPI_ROOT_PROCESS)
        return;
#endif
    SFile *sf = SFile::Create(this->output, SFILE_MODE_WRITE);

    //slibreal_t pixls[3] = {(slibreal_t)pixelsX, (slibreal_t)pixelsY, (slibreal_t)pixelsZ};
    //sf->WriteList("pixels", pixls, 3);
    sf->WriteScalar("pixels", pixelsX);

    //sf->WriteList("image", this->s3dimage, this->imagesize);
    sfilesize_t ndims = 3, dims[3] = {(sfilesize_t)pixelsX, (sfilesize_t)pixelsY, (sfilesize_t)pixelsZ};
    sf->WriteMultiArray("image", this->s3dimage, ndims, dims);

    sf->WriteList("xmin", &(this->point0[0]), 1);
    sf->WriteList("xmax", &(this->point1[0]), 1);
    sf->WriteList("ymin", &(this->point0[1]), 1);
    sf->WriteList("ymax", &(this->point1[1]), 1);
    sf->WriteList("zmin", &(this->point0[2]), 1);
    sf->WriteList("zmax", &(this->point1[2]), 1);

	this->WriteCommonQuantities(sf);

    sf->Close();

    SOFT::PrintInfo("Wrote S3D to '%s'.", this->output.c_str());

    delete [] s3dimage;
}

