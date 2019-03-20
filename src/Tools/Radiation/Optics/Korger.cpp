/**
 * Optics for a simple polarimeter consisting of a linear
 * polarization filter and a quarter-wave plate. Based on
 * the model by
 *
 *   Korger et al., Opt. Express 21, 27032-27042 (2013)
 *   https://doi.org/10.1364/OE.21.027032
 */

#include <complex>
#include <iostream>
#include <softlib/config.h>
#include <softlib/Matrix.h>
#include "Tools/Radiation/Optics/Korger.h"

using namespace std;
using namespace __Radiation;

typedef Matrix<3,3,complex<slibreal_t>> cMatrix;
typedef Vector<3,complex<slibreal_t>> cVector;

/**
 * Constructor. 
 */
Korger::Korger(Detector *det, ConfigBlock *conf) : Optics(det) {
    Configure(conf);
}

/**
 * Configure the Korger model.
 *
 * conf: Configuration block containing settings.
 */
void Korger::Configure(ConfigBlock*) {
    // No options
}

/**
 * Apply the Korger model to a spectrum
 * of electric field components coming
 * from the same source direction.
 */
void Korger::ApplyOptics(
    const struct Optics::Efield &E,
    slibreal_t *I, slibreal_t *Q,
    slibreal_t *U, slibreal_t *V
) {
    // Construct projection operator
    cVector
        ap00 = detector->GetEHat1(),
        ap45 = 1.0/sqrt(2.0) * (detector->GetEHat1() + detector->GetEHat2()),
        ap90 = detector->GetEHat2();
    complex<slibreal_t>
        kap00 = E.zhat.Dot(ap00),
        kap45 = E.zhat.Dot(ap45),
        kap90 = E.zhat.Dot(ap90);
    cVector
        a00 = (ap00 - E.zhat*kap00) / (complex<slibreal_t>)sqrt(1.0 - norm(kap00)),
        a45 = (ap45 - E.zhat*kap45) / (complex<slibreal_t>)sqrt(1.0 - norm(kap45)),
        a90 = (ap90 - E.zhat*kap90) / (complex<slibreal_t>)sqrt(1.0 - norm(kap90));

    cMatrix ONE, Tl;

    // Unit matrix
    ONE(0,0) = ONE(1,1) = ONE(2,2) = 1.0;
    // Quarter-wave retarder
    Tl(0,0) = Tl(2,2) = 1.0;
    Tl(1,1) = -1i;

    cMatrix Tp00 = ONE - cMatrix(a00, a00);
    cMatrix Tp45 = ONE - cMatrix(a45, a45);
    cMatrix Tp90 = ONE - cMatrix(a90, a90);

    /*cVector nhat = E.yhat;
    cVector tv1 = Tp00 * nhat;
    cVector tv2 = Tp90 * nhat;

    slibreal_t x = sqrt(real(tv1.Norm()*tv1.Norm() + tv2.Norm()*tv2.Norm()));
    if (abs(x - 1.0) > 1e-2)
        printf("Wait now... %e !!!\n", abs(x-1.0));*/

    // Construct Stokes parameters from PI(psi)
    for (unsigned int i = 0; i < E.nE; i++) {
        cVector e = E.Ex[i] * E.xhat + E.Ey[i] * E.yhat;

        cVector PI00vec = Tp00 * e;
        cVector PI45vec = Tp45 * e;
        cVector PI90vec = Tp90 * e;
        cVector PIi45vec = (Tp45 * Tl) * e;

        slibreal_t PI00 = real(PI00vec.Norm()); PI00 *= PI00;
        slibreal_t PI45 = real(PI45vec.Norm()); PI45 *= PI45;
        slibreal_t PI90 = real(PI90vec.Norm()); PI90 *= PI90;
        slibreal_t PIi45 = real(PIi45vec.Norm()); PIi45 *= PIi45;

        I[i] = PI00 + PI90;

        if (Q != nullptr) Q[i] = PI00 - PI90;
        if (U != nullptr) U[i] = 2.0 * PI45 - PI00 - PI90;
        if (V != nullptr) V[i] = 2.0 * PIi45 - PI00 - PI90;
    }
}

