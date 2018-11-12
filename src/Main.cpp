/**
 * SOFT
 *
 * Written by: Mathias Hoppe (2018)
 * Email: hoppe@chalmers.se
 */

#include <softlib/Configuration.h>

#include "Init/Init.h"
#include "SOFT.h"
#include "softexit.h"

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
	Configuration *conf;
    SOFT *soft;

	
	/* Handle command-line arguments */
    try {
        conf = new Configuration();
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

	return 0;
}

