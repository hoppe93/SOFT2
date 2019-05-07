/**
 * Optics for a simple polarimeter consisting of a linear
 * polarization filter and a quarter-wave plate. Based on
 * the model by
 *
 *   Korger et al., Opt. Express 21, 27032-27042 (2013)
 *   https://doi.org/10.1364/OE.21.027032
 */

#include <iostream>
#include <softlib/config.h>
#include "Tools/Radiation/Optics/Korger.h"

using namespace std;
using namespace __Radiation;

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
    Vector<3>
        e1 = detector->GetEHat1(),
        e2 = detector->GetEHat2(),
        nh = detector->GetDirection();

    slibreal_t
        a1 = e1.Dot(E.yhat), a2 = e2.Dot(E.yhat),
        b1 = e1.Dot(E.xhat), b2 = e2.Dot(E.xhat);

    for (unsigned int i = 0; i < E.nE; i++) {
        slibreal_t
            PI0  = (E.Ex2[i]*a2*a2 + E.Ey2[i]*b2*b2) / (a2*a2 + b2*b2),
            PI4  = (E.Ex2[i]*(a2-a1)*(a2-a1) + E.Ey2[i]*(b2-b1)*(b2-b1)) / ((a1-a2)*(a1-a2) + (b1-b2)*(b1-b2)),
            PI2  = (E.Ex2[i]*a1*a1 + E.Ey2[i]*b1*b1) / (a1*a1 + b1*b1),
            PIi4 = 0.5*(E.Ex2[i] + E.Ey2[i]) - E.ExEy[i];

        I[i] = PI0 + PI2;

        if (Q != nullptr) Q[i] = PI0 - PI2;
        if (U != nullptr) U[i] = 2.0*PI4  - PI0 - PI2;
        if (V != nullptr) V[i] = 2.0*PIi4 - PI0 - PI2;
    }
}

