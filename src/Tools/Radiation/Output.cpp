/**
 * Implementation of the output routine for the
 * 'Radiation' tool.
 */

#include <softlib/config.h>
#include "Tools/Radiation/Radiation.h"

using namespace __Radiation;

/**
 * Finish the simulation on this thread.
 */
void Radiation::Finish() {
    for (unsigned int i = 0; i < this->noutput; i++) {
        this->output[i]->Finish();
    }
}

/**
 * Generate output.
 */
void Radiation::Output() {
    for (unsigned int i = 0; i < this->noutput; i++) {
        this->output[i]->Generate();
    }
}
