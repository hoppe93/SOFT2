/**
 * Simplest "quadrature" rule implemented. Just evaluates
 * the emission function at the center point of the detector.
 */

#include <softlib/config.h>
#include "Tools/Radiation/Models/AngularDistribution.h"
#include "Tools/Radiation/Models/AngularDistribution/Quadrature2D/ADEval2D.h"

using namespace std;
using namespace __Radiation;

/**
 * Evaluate the function at the center point of the detector.
 */
slibreal_t ADEval2D::Integrate(
    RadiationParticle *rp, bool pol,
    slibreal_t *I, slibreal_t *Q, slibreal_t *U, slibreal_t *V
) {
    struct angles a;

    EvaluateAngles(0.0, 0.0, rp, a);
    
    this->emission->Evaluate(a.n, a.sinMu, a.cosMu, pol);
    slibreal_t area = this->detector->GetAperture();
    area *= area;
    slibreal_t nDotNHat = rp->GetRCPHat().Dot(this->detector->GetDirection());
    slibreal_t fac = area * nDotNHat / a.rcp2;

    slibreal_t pwr = fac * this->emission->GetTotalEmission();

    if (nwavelengths > 0) {
        if (pol)
            MultiplySpectra<true>(
                fac, nwavelengths,
                I, Q, U, V,
                this->emission->GetStokesI(),
                this->emission->GetStokesQ(),
                this->emission->GetStokesU(),
                this->emission->GetStokesV()
            );
        else
            MultiplySpectra<false>(
                fac, nwavelengths,
                I, Q, U, V,
                this->emission->GetStokesI(),
                this->emission->GetStokesQ(),
                this->emission->GetStokesU(),
                this->emission->GetStokesV()
            );
    }

    return pwr;
}

