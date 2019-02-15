/**
 * Optics for a simple polarimeter consisting of a linear
 * polarization filter and a quarter-wave plate. Based on
 * the model by
 *
 *   Korger et al., Opt. Express 21, 27032-27042 (2013)
 *   https://doi.org/10.1364/OE.21.027032
 */

#include <complex>
#include "Tools/Radiation/Optics/Korger.h"

using namespace std;


Korger::Korger(ConfigBlock *conf) {
    Configure(conf);
}

/**
 * Configure the Korger model.
 *
 * conf: Configuration block containing settings.
 */
void Korger::Configure(ConfigBlock*) {
    // No options
}

/**
 * Apply the Korger model to a spectrum
 * of electric field components coming
 * from the same source direction.
 */
void Korger::ApplyOptics(
    const struct Optics::Efield &E,
    slibreal_t *I, slibreal_t *Q,
    slibreal_t *U, slibreal_t *V
) {
    
}

