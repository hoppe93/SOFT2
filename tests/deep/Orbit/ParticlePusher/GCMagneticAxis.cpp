/**
 * Test exact solution for guiding-center equations
 * of motion to zeroth order on magnetic axis.
 */

#include <string>
#include <softlib/config.h>
#include <softlib/constants.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include "Init/InitConfig.h"
#include "PhaseSpace/Particle.h"
#include "SOFT.h"
#include "Test_ParticlePusher.h"

/**
 * Solve the guiding-center equations of motion
 * to zeroth-order for a guiding-center starting on
 * the magnetic axis, and compare to the exact
 * solution.
 *
 * RETURNS true if the ParticlePusher yields the correct
 * solution, false otherwise.
 */
bool Test_ParticlePusher::GCMagneticAxis() {
    slibreal_t
        p = 10.0,
        thetap = 0.1,
        gamma = sqrt(p*p + 1.0),
        ppar = p*cos(thetap);

    Orbit *o = this->GenerateOrbit(Rm, p, thetap, false);

    // Compare to exact solution
    unsigned int ntau = o->GetNTau(), i;
    slibreal_t tau, *X,
               Xex[3], wT = LIGHTSPEED*ppar / (gamma*Rm);

    for (i = 0; i < ntau; i++) {
        tau = o->GetTau(i);
        Xex[0] = Rm * cos(wT*tau);
        Xex[1] =-Rm * sin(wT*tau);
        Xex[2] = 0.0;

        X = o->GetX(i);

        for (int j = 0; j < 2; j++)
            if (fabs((X[j]-Xex[j])/Rm) > Test_ParticlePusher::GCMA_TOLERANCE)
                return false;
    }

    delete o;

    return true;
}

