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
Topview::Topview(
	Detector *d, MagneticField2D *m,
	ParticleGenerator *pgen, SOFT *soft
) : RadiationOutput(d, m, pgen, soft) {
    this->max_radius = m->GetMaxRadius();
}

/**
 * Destructor.
 */
Topview::~Topview() { }

