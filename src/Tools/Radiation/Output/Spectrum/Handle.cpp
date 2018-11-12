/**
 * Radiation :: Output :: Spectrum
 * 
 * Handle output from a RadiationModel.
 */

#include <softlib/config.h>
#include "Tools/Radiation/Output/Spectrum.h"

using namespace __Radiation;

/**
 * Handle output from the given model.
 *
 * det:   Detector observing the radiation.
 * model: Radiation model to handle output from.
 * rp:    Particle producing the radiation.
 */
void Spectrum::Handle(Detector *__UNUSED__(det), Model *model, RadiationParticle *rp) {
    slibreal_t diffel, f;
    slibreal_t *_I, *_Q, *_U, *_V;

    f = rp->GetF();
    diffel = rp->GetDifferentialElement();
    _I = model->GetStokesI();

    for (unsigned int i = 0; i < this->nwavelengths; i++)
        this->I[i] += _I[i] * diffel * f;

    if (this->MeasuresPolarization()) {
        _Q = model->GetStokesQ();
        _U = model->GetStokesU();
        _V = model->GetStokesV();

        for (unsigned int i = 0; i < this->nwavelengths; i++) {
            this->Q[i] += _Q[i] * diffel * f;
            this->U[i] += _U[i] * diffel * f;
            this->V[i] += _V[i] * diffel * f;
        }
    }
}

