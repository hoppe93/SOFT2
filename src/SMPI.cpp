/**
 * Custom MPI wrapper for SOFT
 */

#ifndef WITH_MPI
#   error "The SMPI module requires SOFT to be compiled with MPI."
#endif

#include <mpi.h>
#include "smpi.h"

namespace SMPI {
    MPI_Datatype MPI_SLIBREAL_T;

    /**
     * Initialize MPI for SOFT by making a call to
     * MPI_Init() and register relevant data types.
     *
     * argc: Number of command-line arguments passed to the program.
     * argv: Command-line arguments passed to the program.
     */
    void init(int argc, char *argv[]) {
        MPI_Init(&argc, &argv);

        // Register data types
        if (std::is_same<slibreal_t, float>::value)
            verify(MPI_Type_contiguous(1, MPI_FLOAT, &MPI_SLIBREAL_T));
        else if (std::is_same<slibreal_t, double>::value)
            verify(MPI_Type_contiguous(1, MPI_DOUBLE, &MPI_SLIBREAL_T));
        else if (std::is_same<slibreal_t, long double>::value)
            verify(MPI_Type_contiguous(1, MPI_LONG_DOUBLE, &MPI_SLIBREAL_T));
        else
            throw MPIException("Unable to register MPI datatype for 'slibreal_t': 'slibreal_t' is of an unknown type.");

        verify(MPI_Type_commit(&MPI_SLIBREAL_T));
    }

    void finalize() {
        MPI_Finalize();
    }

    void verify(int err) {
        if (MPI_SUCCESS)
            return;
        else
            throw MPIException(err);
    }
}
