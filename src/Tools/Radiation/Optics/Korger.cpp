/**
 * Optics for a simple polarimeter consisting of a linear
 * polarization filter and a quarter-wave plate. Based on
 * the model by
 *
 *   Korger et al., Opt. Express 21, 27032-27042 (2013)
 *   https://doi.org/10.1364/OE.21.027032
 */

#include <complex>
#include "Tools/Radiation/Optics/Korger.h"

using namespace std;


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
    Vector<3>
        ap00 = detector->ehat1,
        ap45 = 1.0/sqrt(2.0) * (detector->ehat1 + detector->ehat2);
        ap90 = detector->ehat2;
    slibreal_t
        kap00 = E.zhat.Dot(ap00),
        kap45 = E.zhat.Dot(ap45),
        kap90 = E.zhat.Dot(ap90);
    Vector<3,complex<slibreal_t>>
        a00 = (ap00 - E.zhat*kap00) / sqrt(1.0 - kap00*kap00),
        a45 = (ap45 - E.zhat*kap45) / sqrt(1.0 - kap45*kap45),
        a90 = (ap90 - E.zhat*kap90) / sqrt(1.0 - kap90*kap90);

    Matrix<3,3,complex<slibreal_t>>
        ONE, Tl, Tp00, Tp45, Tp90;

    // Unit matrix
    ONE[0,0] = ONE[1,1] = ONE[2,2] = 1.0;
    // Quarter-wave retarder
    Tl[0,0] = Tl[2,2] = 1.0;
    Tl[1,1] = -1i;

    Tp00 = ONE - Matrix<3,3,complex<slibreal_t>>(a00, a00);
    Tp45 = ONE - Matrix<3,3,complex<slibreal_t>>(a45, a45);
    Tp90 = ONE - Matrix<3,3,complex<slibreal_t>>(a90, a90);

    // Create irradiances (PI(psi))
    // Construct Stokes parameters from PI(psi)
    for (unsigned int i = 0; i < nE; i++) {
        Vector<3,complex<slibreal_t>> e = E.Ex[i] * E.xhat + E.Ey[i] + E.yhat;

        Vector<3,complex<slibreal_t>> PI00vec = Tp00 * e;
        Vector<3,complex<slibreal_t>> PI45vec = Tp45 * e;
        Vector<3,complex<slibreal_t>> PI90vec = Tp90 * e;
        Vector<3,complex<slibreal_t>> PIi45vec = (Tp45 * Tl) * e;

        slibreal_t PI00 = real(PI00vec.Norm()); PI00 *= PI00;
        slibreal_t PI45 = real(PI45vec.Norm()); PI45 *= PI45;
        slibreal_t PI90 = real(PI90vec.Norm()); PI90 *= PI90;
        slibreal_t PIi45 = real(PIi45vec.Norm()); PIi45 *= PIi45;

        I[i] += PI00 + PI90;
        Q[i] += PI00 - PI90;
        U[i] += 2.0 * PI45 - PI00 - PI90;
        V[i] += 2.0 * PIi45 - PI00 - PI90;
    }
}

