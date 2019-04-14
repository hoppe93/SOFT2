/**
 * Implementation of the topview output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Output/Topview.h"

using namespace __Radiation;

/**
 * Constructor.
 */
Topview::Topview(Detector *d, MagneticField2D *m, ParticleGenerator *pgen) : RadiationOutput(d, m, pgen) {
    this->max_radius = m->GetMaxRadius();
}

/**
 * Destructor.
 */
Topview::~Topview() { }

