/**
 * Routines for neately exiting SOFT.
 */

#include <cstdlib>

/**
 * Stops execution and exits safely.
 * The main purpose of this function is
 * to allow MPI-safe exits.
 *
 * exitcode: Return code to exit with.
 */
void softexit(int exitcode) {
	/* TODO Handle MPI-safe exit */
	exit(exitcode);
}
