/**
 * Test the 2D quadrature rules used in the angular
 * distribution emission model.
 *
 *   [ ] - Verify that the quadrature rule recovers the
 *         known result for an isotropically emitting
 *         particle.
 */

#include <algorithm>
#include <softlib/config.h>
#include "Tools/Radiation/Detector.h"
#include "Test_ADQuadrature2D.h"

#include "Tools/Radiation/Models/Isotropic.h"
#include "Tools/Radiation/Models/AngularDistribution/Quadrature2D/ADEval2D.h"
#include "Tools/Radiation/Models/AngularDistribution/Quadrature2D/ADSimpson2D.h"

using namespace std;
using namespace __Radiation;

const slibreal_t Test_ADQuadrature2D::DETECTOR_APERTURE=0.006;
const unsigned int Test_ADQuadrature2D::NTESTPARTICLES=2;
const slibreal_t Test_ADQuadrature2D::test_particles[NTESTPARTICLES][3] = {
    {0.68,0.0,0.0},
    {0.75,0.20,0.03}
};

/**
 * Returns a sample detector.
 */
Detector *Test_ADQuadrature2D::GetDetector(unsigned int nwavelengths) {
    const slibreal_t
        aperture = DETECTOR_APERTURE,
        visang   = 1.0,
        dir[3] = {0.0,1.0,0.0},
        pos[3] = {0.0,-1.069,0.0};
        
    Vector<3>
        direction(dir),
        position(pos);

    return new Detector(aperture, visang, direction, position, nwavelengths, 0.1, 0.9);
}

/**
 * Returns the pre-defined RadiationParticle object
 * with index 'i'.
 *
 * i: Index in table 'predef_particles' with pre-defined
 *    particle parameters to test.
 */
RadiationParticle *Test_ADQuadrature2D::GetRadiationParticle(unsigned int i, Detector *det) {
    slibreal_t
        Jdtdrho = 1.0,
        ppar = 1.0, pperp = 1.0,
        gamma = 1.0, p2 = 1.0,
        B = 1.0, m = 1.0, q = 1.0,
        BB[3] = {1.0,0.0,0.0},
        Pp[3] = {1.0,0.0,0.0};
    Vector<3> P(Pp), Bvec(BB), bHat = Bvec / B;

    if (i >= NTESTPARTICLES)
        throw SOFTException("Trying to access non-existant test-particle.");

    return new RadiationParticle(
        test_particles[i], P,
        Jdtdrho, ppar, pperp,
        gamma, p2, det->GetPosition(),
        B, Bvec, bHat, m, q, 0, 0, 0
    );
}

/**
 * Conduct the 'isotropic emission' test.
 */
template<class T, typename ... Args>
bool Test_ADQuadrature2D::IsotropicEmission(slibreal_t tol, Args&& ... args) {
    bool success = true;

    // Test with no spectrum
    Detector *det1 = GetDetector(0);
    Test_ADQuadrature2D_Emission em1(det1);
    T q1(&em1, det1, forward<Args>(args) ...);

    for (unsigned int i = 0; i < NTESTPARTICLES; i++)
        success &= IsotropicEmissionTest(q1, GetRadiationParticle(i, det1), tol);

    // Test with spectrum
    Detector *det2 = GetDetector(40);
    Test_ADQuadrature2D_Emission em2(det2);
    T q2(&em2, det2, forward<Args>(args) ...);

    for (unsigned int i = 0; i < NTESTPARTICLES; i++)
        success &= IsotropicEmissionTest(q2, GetRadiationParticle(i, det2), tol);

    delete det1;
    delete det2;

    return success;
}

/**
 * Test a specific 2D quadrature rule with
 * a specific test particle.
 *
 * q:  2D quadrature rule to test.
 * rp: RadiationParticle object specifying the particle emitting state.
 */
bool Test_ADQuadrature2D::IsotropicEmissionTest(ADQuadrature2D &q, RadiationParticle *rp, slibreal_t tol) {
    bool success = true;
    Detector *det = q.GetDetector();
    slibreal_t
        *I = new slibreal_t[det->GetNWavelengths()],
        *Q = new slibreal_t[det->GetNWavelengths()],
        *U = new slibreal_t[det->GetNWavelengths()],
        *V = new slibreal_t[det->GetNWavelengths()];

    // Initialize
    unsigned int i;
    for (i = 0; i < det->GetNWavelengths(); i++)
        I[i] = Q[i] = U[i] = V[i] = 0.0;

    // Check...
    success &= IsotropicEmissionVerify(q, rp, false, tol, I, Q, U, V);
    if (det->GetNWavelengths() > 0)
        success &= IsotropicEmissionVerify(q, rp, true,  tol, I, Q, U, V);

    delete [] I;
    delete [] Q;
    delete [] U;
    delete [] V;

    return success;
}

/**
 * Verify that the given quadrature rule returns
 * the correct result for the given RadiationParticle.
 */
bool Test_ADQuadrature2D::IsotropicEmissionVerify(
    ADQuadrature2D &q, RadiationParticle *rp,
    bool withPolarization, slibreal_t tol,
    slibreal_t *I, slibreal_t *Q, slibreal_t *U, slibreal_t *V
) {
    // Calculate expected result
    Detector *det = q.GetDetector();
    Vector<3>
        &nhat = det->GetDirection(),
        &e1 = det->GetEHat1(),
        &e2 = det->GetEHat2(),
        &rcp = rp->GetRCP();

    const slibreal_t
        cn = nhat.Dot(rcp),
        c1 = e1.Dot(rcp),
        c2 = e2.Dot(rcp),
        l  = 0.5 * det->GetAperture(),
        correct =
            Isotropic::LambdaIso( l,  l, cn, c1, c2) +
            Isotropic::LambdaIso(-l, -l, cn, c1, c2) -
            Isotropic::LambdaIso(-l,  l, cn, c1, c2) -
            Isotropic::LambdaIso( l, -l, cn, c1, c2);

    slibreal_t pwr = q.Integrate(rp, withPolarization, I, Q, U, V);

    // First, check returned power...
    slibreal_t Delta = fabs((pwr-correct)/correct);
    if (isnan(Delta) || Delta >= tol) {
        this->PrintError("Integrated (total) isotropic emission does not match with theory. Delta = %e", Delta);
        return false;
    }

    // Next, check spectrum and polarization components
    // if calculated.
    if (det->GetNWavelengths() > 0) {
        unsigned int i, n = det->GetNWavelengths();
        slibreal_t c;
        for (i = 0; i < n; i++) {
            c = (i+1)*correct;
            Delta = fabs((I[i]-c)/c);

            if (withPolarization) {
                Delta = max(Delta, fabs((Q[i]-c)/c));
                Delta = max(Delta, fabs((U[i]-c)/c));
                Delta = max(Delta, fabs((V[i]-c)/c));
            }

            if (isnan(Delta) || Delta >= tol) {
                this->PrintError("Integrated (spectral) emission does not match at wavelength #%u. Max Delta = %e", Delta, i+1);
                return false;
            }
        }
    }

    return true;
}

/**
 * Run all tests of this module.
 */
bool Test_ADQuadrature2D::Run(bool) {
    bool success = true;
    slibreal_t nSimpsonPoints = 10;
    slibreal_t SIMPSON2D_TOL = 1e-5;

    if (!IsotropicEmission<ADEval2D>(SIMPSON2D_TOL))
        success = false;
    if (!IsotropicEmission<ADSimpson2D>(SIMPSON2D_TOL, nSimpsonPoints))
        success = false;

    return success;
}

