/**
 * Routines for neately exiting SOFT.
 */

#include <cstdlib>
#include "config.h"

#ifdef WITH_MPI
#   include <mpi.h>
#   include "SMPI.h"
#endif

/**
 * Stops execution and exits safely.
 * The main purpose of this function is
 * to allow MPI-safe exits.
 *
 * exitcode: Return code to exit with.
 */
void softexit(int exitcode) {
#ifdef WITH_MPI
    SMPI::finalize();
#endif

	exit(exitcode);
}
