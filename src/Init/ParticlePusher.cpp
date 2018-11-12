/**
 * Handles default settings for the 'ParticlePusher'
 * class.
 */

#include <string>
#include "Orbit/ParticlePusher.h"

using namespace std;

/**
 * Constructs and returns a SOFT configuration string
 * that configures global defaults for the ParticlePusher.
 *
 * RETURNS a string that can be interpreted by a
 * 'Configuration' object.
 */
string InitPusherDefaults() {
    return ParticlePusher::equation_defaults;
}
