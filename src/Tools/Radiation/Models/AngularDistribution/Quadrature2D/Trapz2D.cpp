/**
 * Radiation model 'AngularDistribution'
 *
 * Implementation of a 2D trapezoidal quadrature rule
 * for the integration over the detector surface.
 */

#include <softlib/config.h>
#include "Tools/Radiation/Models/AngularDistribution.h"
#include "Tools/Radiation/Models/AngularDistributionException.h"
#include "Tools/Radiation/Models/AngularDistribution/Quadrature2D/ADTrapz2D.h"

using namespace __Radiation;

/**
 * Constructor.
 */
ADTrapz2D::ADTrapz2D(ADEmission *em, Detector *det, unsigned int nsamples) : ADQuadrature2D(em, det) {
    this->nsamples = nsamples;
}

/**
 * Carry out a 2D trapezoidal integration over the detector surface.
 *
 * rp: Object representing the particle emission state.
 */
slibreal_t ADTrapz2D::Integrate(
    RadiationParticle *rp, bool pol,
    slibreal_t *I, slibreal_t *Q, slibreal_t *U, slibreal_t *V
) {
    if (nwavelengths == 0)
        return _OuterIntegral<false, false>(rp, I, Q, U, V);
    else {
        if (!pol)
            return _OuterIntegral<true, false>(rp, I, Q, U, V);
        else
            return _OuterIntegral<true, true>(rp, I, Q, U, V);
    }
}

/**
 * Outer integral.
 *
 * rp:         Particle emitting state.
 * I, Q, U, V: Arrays to store the integrated Stokes parameters in.
 */
template<bool withSpectrum, bool withPolarization>
slibreal_t ADTrapz2D::_OuterIntegral(
    RadiationParticle *rp,
    slibreal_t *I, slibreal_t *Q, slibreal_t *U, slibreal_t *V
) {
    slibreal_t S, dX,
        rdet  = this->detector->GetAperture(),
        rdet2 = 0.5 * rdet,
        // The detector direction is orthogonal to e1 & e2, so we don't
        // need to reevaluate this product for each new rcp in the integral
        nDotNHat = this->detector->GetDirection().Dot(rp->GetRCPHat());
    unsigned int i;

    dX = rdet / (nsamples-1);

    if (withSpectrum)
        ResetSpectra<withPolarization>(nwavelengths, this->I, this->Q, this->U, this->V);

    // First endpoint
    S  = 0.5*_InnerIntegral<withSpectrum, withPolarization>(rdet2, dX, rdet2, rp);
    if (withSpectrum)
        SumSpectra<withPolarization>(
            0.5, nwavelengths, I, Q, U, V,
            this->I, this->Q, this->U, this->V
        );

    // Second endpoint
    S += 0.5*_InnerIntegral<withSpectrum, withPolarization>(-rdet2, dX, rdet2, rp);
    if (withSpectrum)
        SumSpectra<withPolarization>(
            0.5, nwavelengths, I, Q, U, V,
            this->I, this->Q, this->U, this->V
        );

    // Inner points
    for (i = 1; i < nsamples-1; i++) {
        S += _InnerIntegral<withSpectrum, withPolarization>(i*dX - rdet2, dX, rdet2, rp);
        if (withSpectrum)
            SumSpectra<withPolarization>(
                1.0, nwavelengths, I, Q, U, V,
                this->I, this->Q, this->U, this->V
            );
    }

    // Multiply with the differential element
    slibreal_t d = nDotNHat * dX*dX;
    if (withSpectrum)
        MultiplySpectra<withPolarization>(
            d, nwavelengths, I, Q, U, V
        );

    return (S * d);
}

/**
 * Inner integral.
 *
 * X: X-coordinate of integration.
 */
#define EVAL(s,p) \
    (s? \
        (p? (this->emission->CalculatePolarization(a.n, a.sinMu, a.cosMu)):\
            (this->emission->CalculateSpectrum(a.n, a.sinMu, a.cosMu))\
        ):\
    this->emission->CalculateAngularDistribution(a.n, a.sinMu, a.cosMu))

template<bool withSpectrum, bool withPolarization>
slibreal_t ADTrapz2D::_InnerIntegral(
    slibreal_t X, slibreal_t dX, slibreal_t rdet2,
    RadiationParticle *rp
) {
    slibreal_t S;
    struct angles a;
    unsigned int i;

    if (withSpectrum)
        ResetSpectra<withPolarization>(
            nwavelengths, I, Q, U, V
        );

    // Endpoint 1
    EvaluateAngles(X, rdet2, rp, a);
    EVAL(withSpectrum, withPolarization);
    S  = 0.5*this->emission->GetTotalEmission() / a.rcp2;
    if (withSpectrum)
        SumSpectra<withPolarization>(
            0.5 / a.rcp2, nwavelengths, I, Q, U, V,
            this->emission->GetStokesI(),
            this->emission->GetStokesQ(),
            this->emission->GetStokesU(),
            this->emission->GetStokesV()
        );

    // Endpoint 2
    EvaluateAngles(X, -rdet2, rp, a);
    EVAL(withSpectrum, withPolarization);
    S += 0.5*this->emission->GetTotalEmission() / a.rcp2;
    if (withSpectrum)
        SumSpectra<withPolarization>(
            0.5 / a.rcp2, nwavelengths, I, Q, U, V,
            this->emission->GetStokesI(),
            this->emission->GetStokesQ(),
            this->emission->GetStokesU(),
            this->emission->GetStokesV()
        );

    // Inner points
    for (i = 1; i < nsamples-1; i++) {
        EvaluateAngles(X, i*dX - rdet2, rp, a);
        EVAL(withSpectrum, withPolarization);

        S += this->emission->GetTotalEmission() / a.rcp2;
        if (withSpectrum)
            SumSpectra<withPolarization>(
                1.0 / a.rcp2, nwavelengths, I, Q, U, V,
                this->emission->GetStokesI(),
                this->emission->GetStokesQ(),
                this->emission->GetStokesU(),
                this->emission->GetStokesV()
            );
    }

    return S;
}

