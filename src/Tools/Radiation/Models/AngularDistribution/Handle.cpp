/**
 * Handle a single point of phase-space with the angular distribution model.
 */

#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/AngularDistribution.h"

using namespace __Radiation;

/**
 * Constructor.
 */
AngularDistribution::AngularDistribution(Radiation *rad) : Model(rad) {
    unsigned int n = rad->detector->GetNWavelengths(), i;

    if (n > 0) {
        this->I = new slibreal_t[n];
        this->Q = new slibreal_t[n];
        this->U = new slibreal_t[n];
        this->V = new slibreal_t[n];

        for (i = 0; i < n; i++)
            I[i] = Q[i] = U[i] = V[i] = 0.0;
    } else
        this->I = this->Q = this->U = this->V = nullptr;
}
AngularDistribution::~AngularDistribution() {
    delete emission;

    if (this->I != nullptr)
        delete [] this->I;
    if (this->Q != nullptr)
        delete [] this->Q;
    if (this->U != nullptr)
        delete [] this->U;
    if (this->V != nullptr)
        delete [] this->V;
}

/**
 * Calculates the amount of emitted radiation in this
 * timestep.
 *
 * rp: RadiationParticle defining the instantaneous
 *     particle properties.
 */
void AngularDistribution::InitializeTimestep(RadiationParticle *rp) {
    this->emission->Prepare(rp, this->parent->MeasuresPolarization());
}

/**
 * Handles the given particle radiation state and
 * calculates the amount of radiation reaching
 * the detector.
 *
 * rp:             RadiationParticle defining the instantaneous
 *                 particle properties.
 * sinphi, cosphi: Sine/cosine of current toroidal angle.
 */
void AngularDistribution::HandleParticle(RadiationParticle *rp, const slibreal_t sinphi, const slibreal_t cosphi) {
    this->emission->InitializeToroidalStep(sinphi, cosphi);

    this->quadrature2d->Integrate(
        rp, this->parent->MeasuresPolarization(),
        this->I, this->Q, this->U, this->V
    );

    this->nonzero = (this->emission->GetTotalEmission() > 0);
}

