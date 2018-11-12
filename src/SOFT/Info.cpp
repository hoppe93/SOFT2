/**
 * Module for printing information about a SOFT
 * object.
 */

#include <string>
#include <vector>
#include "config.h"
#include "SOFTLocal.h"

const char soft_git_refspec[] = SOFT_GIT_REFSPEC;
const char soft_git_sha1[] = SOFT_GIT_SHA1;

using namespace std;

/**
 * Prints 'welcome' message, showing how the SOFT
 * object has been initialized.
 */
void SOFTLocal::Welcome() {
    struct global_settings *globset = this->soft->GetGlobalSettings();

    // Phase-space size
    unsigned int nr = this->partgen->GetNr(),
                 n1 = this->partgen->GetN1(),
                 n2 = this->partgen->GetN2();

    // Gather tool names
    const vector<Tool*> ltools = thandler->GetTools();
    string tools, t;
    unsigned int ntools = ltools.size(), i;
    for (i = 0; i < ntools; i++) {
        t = "@" + ltools[i]->GetType() + "(" + ltools[i]->GetName() + ")";
        if (i == 0)
            tools = t;
        else
            tools += ", " + t;
    }

    // Print info
    SOFT::PrintInfo("This is SOFT v2 (commit %s)\n", soft_git_sha1);

    SOFT::PrintInfo("No. threads:         %u", globset->num_threads);
    SOFT::PrintInfo("Magnetic field:      %s", globset->magnetic_field.c_str());
    SOFT::PrintInfo("Particle generator:  %s", globset->particle_generator.c_str());
    SOFT::PrintInfo("Enabled tools (%2u):  %s", ntools, tools.c_str());

    SOFT::PrintInfo("Size of phase-space:    NR x    N1 x    N2 =");
    SOFT::PrintInfo("                   = %5u x %5u x %5u = %u", nr, n1, n2, nr*n1*n2);

    SOFT::PrintInfo();

    // Print Tool introductions
    for (i = 0; i < ntools; i++) {
#ifdef COLOR_TERMINAL
        SOFT::PrintInfo("\e[1m:: "+ltools[i]->GetTypeAndName()+"\e[0m");
#else
        SOFT::PrintInfo(":: "+ltools[i]->GetTypeAndName());
#endif
        ltools[i]->Welcome();
    }
}

