/**
 * Radiation model 'AngularDistribution'
 *
 * Implementation of a 2D Simpson quadrature rule
 * for the integration over the detector surface.
 */

#include <softlib/config.h>
#include "Tools/Radiation/Models/AngularDistribution.h"
#include "Tools/Radiation/Models/AngularDistributionException.h"
#include "Tools/Radiation/Models/AngularDistribution/Quadrature2D/ADSimpson2D.h"

using namespace __Radiation;

/**
 * Constructor.
 */
ADSimpson2D::ADSimpson2D(ADEmission *em, Detector *det, unsigned int nsamples) : ADQuadrature2D(em, det) {
    if (nsamples % 2 == 1)
        throw AngularDistributionException("The number of points of the quadrature must be even when using the 'simpson' quadrature.");

    this->nsamples = nsamples;
}

/**
 * Carry out a 2D Simpson integration over the detector surface.
 *
 * rp: Object representing the particle emission state.
 */
slibreal_t ADSimpson2D::Integrate(
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
slibreal_t ADSimpson2D::_OuterIntegral(
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

    dX = rdet / nsamples;

    if (withSpectrum)
        ResetSpectra<withPolarization>(nwavelengths, this->I, this->Q, this->U, this->V);

    // First endpoint
    S  = _InnerIntegral<withSpectrum, withPolarization>(rdet2, dX, rdet2, rp);
    if (withSpectrum)
        SumSpectra<withPolarization>(
            1.0, nwavelengths, I, Q, U, V,
            this->I, this->Q, this->U, this->V
        );

    // Second endpoint
    S += _InnerIntegral<withSpectrum, withPolarization>(-rdet2, dX, rdet2, rp);
    if (withSpectrum)
        SumSpectra<withPolarization>(
            1.0, nwavelengths, I, Q, U, V,
            this->I, this->Q, this->U, this->V
        );

    // Inner points
    for (i = 1; i < nsamples; i += 2) {
        S += 4.0 * _InnerIntegral<withSpectrum, withPolarization>(i*dX - rdet2, dX, rdet2, rp);
        if (withSpectrum)
            SumSpectra<withPolarization>(
                4.0, nwavelengths, I, Q, U, V,
                this->I, this->Q, this->U, this->V
            );
    }
    for (i = 2; i < nsamples-1; i += 2) {
        S += 2.0 * _InnerIntegral<withSpectrum, withPolarization>(i*dX - rdet2, dX, rdet2, rp);
        if (withSpectrum)
            SumSpectra<withPolarization>(
                2.0, nwavelengths, I, Q, U, V,
                this->I, this->Q, this->U, this->V
            );
    }

    // Account for the differential element
    // and factor 1/3 in the Y integral as well
    slibreal_t d = nDotNHat * dX*dX / 9.0;
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
        (p? (this->emission->CalculatePolarization(rp, a.n, a.sinMu, a.cosMu)):\
            (this->emission->CalculateSpectrum(a.n, a.sinMu, a.cosMu))\
        ):\
    this->emission->CalculateAngularDistribution(a.n, a.sinMu, a.cosMu))

template<bool withSpectrum, bool withPolarization>
slibreal_t ADSimpson2D::_InnerIntegral(
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
    //EVAL(withSpectrum, withPolarization);
	this->emission->Evaluate(rp, a.n, a.sinMu, a.cosMu, withPolarization);
    S  = this->emission->GetTotalEmission() / a.rcp2;
    if (withSpectrum)
        SumSpectra<withPolarization>(
            1.0 / a.rcp2, nwavelengths, I, Q, U, V,
            this->emission->GetStokesI(),
            this->emission->GetStokesQ(),
            this->emission->GetStokesU(),
            this->emission->GetStokesV()
        );

    // Endpoint 2
    EvaluateAngles(X, -rdet2, rp, a);
    //EVAL(withSpectrum, withPolarization);
	this->emission->Evaluate(rp, a.n, a.sinMu, a.cosMu, withPolarization);
    S += this->emission->GetTotalEmission() / a.rcp2;
    if (withSpectrum)
        SumSpectra<withPolarization>(
            1.0 / a.rcp2, nwavelengths, I, Q, U, V,
            this->emission->GetStokesI(),
            this->emission->GetStokesQ(),
            this->emission->GetStokesU(),
            this->emission->GetStokesV()
        );

    // Inner points
    for (i = 1; i < nsamples; i += 2) {
        EvaluateAngles(X, i*dX - rdet2, rp, a);
        //EVAL(withSpectrum, withPolarization);
		this->emission->Evaluate(rp, a.n, a.sinMu, a.cosMu, withPolarization);

        S += 4.0 * this->emission->GetTotalEmission() / a.rcp2;
        if (withSpectrum)
            SumSpectra<withPolarization>(
                4.0 / a.rcp2, nwavelengths, I, Q, U, V,
                this->emission->GetStokesI(),
                this->emission->GetStokesQ(),
                this->emission->GetStokesU(),
                this->emission->GetStokesV()
            );
    }

    for (i = 2; i < nsamples-1; i += 2) {
        EvaluateAngles(X, i*dX - rdet2, rp, a);
        //EVAL(withSpectrum, withPolarization);
		this->emission->Evaluate(rp, a.n, a.sinMu, a.cosMu, withPolarization);

        S += 2.0 * this->emission->GetTotalEmission() / a.rcp2;
        if (withSpectrum)
            SumSpectra<withPolarization>(
                2.0 / a.rcp2, nwavelengths, I, Q, U, V,
                this->emission->GetStokesI(),
                this->emission->GetStokesQ(),
                this->emission->GetStokesU(),
                this->emission->GetStokesV()
            );
    }

    return S;
}

