/**
 * Implementation of the unit test for the
 * 'Isotropic' class in SOFT.
 *
 */

#include <cmath>

// From SOFTLib
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>

// From SOFT
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/Models/Isotropic.h"
#include "SOFTException.h"

// From SOFT tests
#include "Test_Isotropic.h"

using namespace std;
using namespace __Radiation;

/**
 * Compare the analytically evaluated isotropic emission
 * function implemented in SOFT to it's numerically
 * evaluated counter-part.
 *
 * iso: Isotropic model to use.
 */
bool Test_Isotropic::CheckIsotropicEmission(Isotropic *iso) {
    slibreal_t x[3] = {0.75,0.20,0.03},
               p[3] = {0.0,-1.0, 0.0 },
               q = -ELEMENTARY_CHARGE,
               m = ELECTRON_MASS,
               B = 1,
               Bvec[3] = {1.0,0.0,0.0};
    Detector *det = iso->GetParent()->detector;
    RadiationParticle *rp = new RadiationParticle(
        x, p, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, det->GetPosition(), B, Bvec, Bvec, m, q, 0, 0, 0
    );

    iso->HandleParticle(rp, 0.0, 1.0);

    slibreal_t
        In = NumericIntegral(det, rp),
        Ia = iso->GetPower(),
        Delta = fabs((In-Ia)/In);

    if (Delta > ISO_TOL) {
        this->PrintError("Isotropic emission does not match the numeric estimate. (Delta = %e)", Delta);
        return false;
    }

    return true;
}

/**
 * Evaluate numerically the isotropic emission integral
 * over the detector surface for the given detector.
 *
 * det: Detector to evaluate the integral for.
 * rp:  Particle to evaluate the integral for.
 */
slibreal_t Test_Isotropic::NumericIntegral(Detector *det, RadiationParticle *rp) {
    Vector<3>
        r0 = rp->GetRCP(),
        nhat  = det->GetDirection(),
        e1 = det->GetEHat1(),
        e2 = det->GetEHat2(),
        X;

    unsigned int i;
    slibreal_t I, h, aperture = det->GetAperture(), _l = 0.5*aperture;
    h = aperture / (slibreal_t)NSIMPSON_STEPS;

    I = NumericIntegral_inner(-_l, _l, h, r0, e1, e2, nhat) +
        NumericIntegral_inner(+_l, _l, h, r0, e1, e2, nhat);

    // Simpson's rule
    for (i = 1; i < NSIMPSON_STEPS; i += 2)
        I += 4.0 * NumericIntegral_inner(i*h-_l, aperture, h, r0, e1, e2, nhat);
    for (i = 2; i < NSIMPSON_STEPS-1; i += 2)
        I += 2.0 * NumericIntegral_inner(i*h-_l, aperture, h, r0, e1, e2, nhat);
    
    return (I * h / 3.0);
}

/**
 * Inner integral in the double integral over the
 * detector surface.
 *
 * l1:   Displacement in ehat1 direction on detector surface.
 * l:    Half detector aperture (integral goes between -l to +l).
 * r0:   rcp relative to central detector position.
 * e1:   Detector plane direction 1.
 * e2:   Detector plane direction 2.
 * nhat: Detector surface normal.
 */
slibreal_t Test_Isotropic::NumericIntegral_inner(
    slibreal_t l1, slibreal_t l, slibreal_t h,
    Vector<3>& r0, Vector<3>& e1, Vector<3>& e2, Vector<3>& nhat
) {
    slibreal_t I = 0;
    unsigned int i;

    I = NumericIntegral_intg(l1, -l, r0, e1, e2, nhat) +
        NumericIntegral_intg(l1, +l, r0, e1, e2, nhat);

    for (i = 1; i < NSIMPSON_STEPS; i += 2)
        I += 4.0 * NumericIntegral_intg(l1, i*h-l, r0, e1, e2, nhat);
    for (i = 2; i < NSIMPSON_STEPS-1; i += 2)
        I += 2.0 * NumericIntegral_intg(l1, i*h-l, r0, e1, e2, nhat);

    return (I * h / 3.0);
}

/**
 * Evaluate the integrand for use in the numeric
 * integral. The integrand is 'n . ^n / r^2'.
 *
 * l1:   Displacement in ehat1 direction on detector surface.
 * l2:   Displacement in ehat2 directoin on detector surface.
 * r0:   rcp relative to central detector position.
 * e1:   Detector plane direction 1.
 * e2:   Detector plane direction 2.
 * nhat: Detector surface normal.
 */
slibreal_t Test_Isotropic::NumericIntegral_intg(
    slibreal_t l1, slibreal_t l2,
    Vector<3>& r0, Vector<3>& e1, Vector<3>& e2, Vector<3>& nhat
) {
    Vector<3> r = r0 - l1*e1 - l2*e2;
    slibreal_t rmag = r.Norm();

    return r.Dot(nhat) / (rmag*rmag*rmag);
}

/**
 * Run the unit test.
 */
bool Test_Isotropic::Run(bool) {
    const slibreal_t detdir[3] = {0.0,1.0,0.0},
               detpos[3] = {0.0,-1.069,0.0};
    Vector<3> vdir(detdir);
    Vector<3> vpos(detpos);

    bool retval = true;
    //MagneticFieldAnalytical2D *mfa = GetMagneticField();

    //Radiation *rad = new Radiation(mfa, nullptr, nullptr);
    Radiation *rad = new Radiation();
    Detector *det = new Detector(
        0.006, 1.1, vdir, vpos, 0
    );
    rad->SetDetector(det);

    Isotropic *iso = new Isotropic(rad);
    iso->SetRadiationValue(1.0);

    if (!CheckIsotropicEmission(iso))
        retval = false;
    else
        this->PrintOK("Isotropic emission function is correct.");

    delete det;

    return retval;
}

