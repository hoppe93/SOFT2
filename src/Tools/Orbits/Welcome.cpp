/**
 * Print a welcome message for the 'Orbits' tool.
 */

#include <string>
#include "SOFT.h"
#include "Tools/Orbits.h"

using namespace std;

void Orbits::Welcome(const string &prefix) {
    SOFT::PrintInfo(prefix+"Not much to say yet...");
}
