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
Spectrum::Spectrum(
	Detector *d, MagneticField2D *m,
	ParticleGenerator *pgen, SOFT *soft
) : RadiationOutput(d, m, pgen, soft) {
    this->nwavelengths = d->GetNWavelengths();

    if (this->nwavelengths == 0)
        throw SpectrumException("Spectrum requested, but detector has no spectral range.");

    this->wavelengths = d->GetWavelengths();
    this->I = new slibreal_t[this->nwavelengths];

    for (unsigned int i = 0; i < this->nwavelengths; i++)
        this->I[i] = 0.0;

    if (this->MeasuresPolarization()) {
        this->Q = new slibreal_t[this->nwavelengths];
        this->U = new slibreal_t[this->nwavelengths];
        this->V = new slibreal_t[this->nwavelengths];

        for (unsigned int i = 0; i < this->nwavelengths; i++) {
            this->Q[i] = 0.0;
            this->U[i] = 0.0;
            this->V[i] = 0.0;
        }
    }
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

