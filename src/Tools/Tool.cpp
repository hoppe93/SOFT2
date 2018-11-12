/**
 * Implementation of a general SOFT tool.
 */

#include <string>
#include "Orbit/Orbit.h"
#include "Tools/Tool.h"

using namespace std;

/**
 * Constructor.
 *
 * name: Name of this tool.
 */
Tool::Tool(const string type) {
    this->type = type;
}

