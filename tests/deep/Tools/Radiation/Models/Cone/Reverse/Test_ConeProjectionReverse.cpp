/**
 * Specific test of reverse cone model.
 */

#include <softlib/Vector.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/Models/Cone/Projection/Reverse.h"
#include "Test_ConeProjectionReverse.h"

using namespace __Radiation;

/**
 * Run the test.
 */
bool Test_ConeProjectionReverse::Run(bool) {
    slibreal_t _nhat[3] = {0.0,1.0,0.0},
               _Rd[3]   = {1.7, -1.7, 0},
               _X[3]    = {1.0005876600749035e+00,1.3633667496475090e+00,0.0000000000000000e+00},
               _phat[3] = {-3.1915546774339076e-01,9.4770237279933911e-01,0.0000000000000000e+00},
               tanThetap = 1.0033467208545054e-01,
               ppar = 10.0,
               pperp = 10.0*tanThetap,
               p2 = ppar*ppar + pperp*pperp,
               gamma = sqrt(1 + p2),
               m = 9.109e-31,
               q = 1.602e-19,
               BB[3] = {1.0,0.0,0.0};
    Vector<3> nhat(_nhat), Rd(_Rd), X(_X), phat(_phat), p;

    Detector *det = new Detector(0.006, 1.0, nhat, Rd, 0);
    ConeProjectionReverse *cp = new ConeProjectionReverse(det);

    RadiationParticle *rp = new RadiationParticle(
        X, phat, 1.0, ppar, pperp, gamma, p2, Rd, 1.0, BB, BB, m, q, 0, 0, 0
    );

    printf("Fraction = %.16e\n", cp->ComputeOverlappingFraction(rp));

    return true;
}
