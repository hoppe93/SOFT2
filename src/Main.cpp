/**
 * SOFT
 *
 * Written by: Mathias Hoppe (2018)
 * Email: hoppe@chalmers.se
 */

#include <softlib/Configuration/ConfigurationScript.h>

#include "config.h"
#include "Init/Init.h"
#include "SOFT.h"
#include "softexit.h"

#ifdef WITH_MPI
#   include <mpi.h>
#   include "SMPI.h"
#endif

#if !defined(NDEBUG) && defined(__linux__)
#	include <fenv.h>
#endif

using namespace std;

/**
 * Program entry point.
 *
 * SOFT takes 0 or 1 input arguments. If one argument is given,
 * it is interpreted as the name of a configuration file to be
 * read. If no argument is given, SOFT reads the configuration
 * file from stdin.
 *
 * Passing more than one argument results in error.
 */
int main(int argc, char *argv[]) {
	ConfigurationScript *conf;
    SOFT *soft;

#if !defined(NDEBUG) && defined(__linux__)
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#ifdef WITH_MPI
    try {
        SMPI::init(&argc, &argv);
    } catch (SMPI::MPIException &ex) {
        SOFT::PrintError("%s", ex.what());
        softexit(1);
    }
#endif
	
	/* Handle command-line arguments */
    try {
        conf = new ConfigurationScript();
        InitConfig(conf);

        if (argc == 1) {
            /* Read config from stdin */
            conf->FromStdin();
        } else if (argc == 2) {
            /* Read config from file */
            string fname = string(argv[1]);
            conf->FromFile(fname);
        } else {
            SOFT::PrintError("Too many command-line arguments.");
            softexit(-1);
        }

        if (conf->HasError())
            softexit(1);

        soft = InitSOFT(conf);
    } catch (SOFTLibException& ex) {
        SOFT::PrintError("%s", ex.what());
        softexit(1);
    }

    try {
        soft->Run();
    } catch (SOFTLibException& ex) {
        SOFT::PrintError("%s", ex.what());
        softexit(1);
    }

#ifdef WITH_MPI
    SMPI::finalize();
#endif

    return 0;
}

