/**
 * Implementation of the spectrum output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Output/Spectrum.h"

using namespace __Radiation;

/**
 * Constructor.
 */
Spectrum::Spectrum(Detector *d, MagneticField2D *m, ParticleGenerator *pgen) : RadiationOutput(d, m, pgen) {
    this->nwavelengths = d->GetNWavelengths();

    if (this->nwavelengths == 0)
        throw SpectrumException("Spectrum requested, but detector has no spectral range.");

    this->wavelengths = d->GetWavelengths();
}

/**
 * Destructor.
 */
Spectrum::~Spectrum() {
    if (this->MeasuresPolarization()) {
        delete [] this->V;
        delete [] this->U;
        delete [] this->Q;
    }

    delete [] this->I;
}

