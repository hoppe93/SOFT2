/**
 * Implementation of the image output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "PhaseSpace/ParticleGenerator.h"
#include "Tools/Radiation/Output/Image.h"

using namespace __Radiation;

/**
 * Constructor.
 */
Image::Image(Detector *d, MagneticField2D *m, ParticleGenerator *__UNUSED__(pgen)) : RadiationOutput(d, m) { }

/**
 * Destructor.
 */
Image::~Image() { }

