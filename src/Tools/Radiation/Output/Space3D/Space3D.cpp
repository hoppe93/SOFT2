/**
 * Implementation of the space3d output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Output/Space3D.h"

using namespace __Radiation;

/**
 * Constructor.
 */
Space3D::Space3D(Detector *d, MagneticField2D *m, ParticleGenerator *pgen) : RadiationOutput(d, m, pgen) { }

/**
 * Destructor.
 */
Space3D::~Space3D() { }

