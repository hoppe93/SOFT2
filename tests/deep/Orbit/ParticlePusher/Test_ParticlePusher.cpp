/**
 * Test the ParticlePusher.
 *
 * Included tests:
 *
 * GUIDING-CENTER
 *   [X] Zeroth-order guiding-center orbit on magnetic axis
 *       is a circle.
 *   [ ] Check that particle energy and magnetic moment are conserved.
 *   [ ] Compare against SOFT v1 guiding-center orbits.
 * PARTICLE
 *   [ ] Comparison to exact solution in straight magnetic field.
 *   [ ] Compare against SOFT v1 particle orbits.
 */

#include <string>
#include <sstream>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/constants.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include "Orbit/ParticlePusher.h"
#include "../../UnitTest.h"
#include "Init/InitConfig.h"
#include "Test_ParticlePusher.h"

/**
 * Generate an Orbit object for a particle with the
 * given momentum and pitch angle.
 * 
 * r:      Particle initial position.
 * p:      Particle momentum.
 * thetap: Particle pitch angle.
 * drifts: If true, includes drifts.
 * eqn:    Name of equation to solve (either
 *         "particle" or "guiding-center").
 */
Orbit *Test_ParticlePusher::GenerateOrbit(const slibreal_t r, const slibreal_t p, const slibreal_t thetap, bool drifts, const string& eqn) {
    ParticlePusher *pp = GenerateOrbit_push(r, p, thetap, drifts, eqn);
    Orbit *o = pp->GetOrbit();

    delete pp;

    return o;
}

ParticlePusher *Test_ParticlePusher::GenerateOrbit_push(const slibreal_t r, const slibreal_t p, const slibreal_t thetap, bool drifts, const string& eqn) {
    slibreal_t
        gamma = sqrt(p*p+1.0),
        ppar = p*cos(thetap);

    // Generate the type of magnetic field implemented
    // in SOFTv1.
    MagneticFieldAnalytical2D *mfa =
        new MagneticFieldAnalytical2D(
            B0, Rm, rminor, MFAFS_CCW, MFAFS_CCW, MFASF_CONSTANT, 1.0, 0.0
        );

    struct global_settings globset;
    globset.include_drifts = drifts;
    globset.particle_pusher = "pusher";

    slibreal_t time = 2*M_PI*r*gamma / (LIGHTSPEED*ppar);
    std::ostringstream ss1, ss2;
    ss1 << std::scientific;
    ss1 << time;

    string ts = ss1.str();
    string config = "equation="+eqn+";\n";
    config += "timeunit=seconds;\n";
    config += "time="+ts+";\n";
    config += "nt=100;\n";

    ss2 << std::scientific;
    ss2 << ORBIT_TOLERANCE;
    
    string tol = ss2.str();
    string eqconfig  = "@Equation guiding-center (guiding-center) {\n";
    eqconfig += "    method = rkdp45;\n";
    eqconfig += "    tolerance = "+tol+";\n";

    if (drifts)
        eqconfig += "    drifts = yes;\n";
    else
        eqconfig += "    drifts = no;\n";

    eqconfig += "}";

    Configuration *conf1 = new Configuration();
    Configuration *conf2 = new Configuration();
    CONFBLOCK_EQUATION = conf2->RegisterBlockType("@EquationGuiding");

    ConfigBlock cb = conf1->FromString(config, "config");
    ConfigBlock eq = conf2->FromString(eqconfig, "eqndefaults");

    ParticlePusher *pp = new ParticlePusher(mfa, &globset, &cb, &eq);

    Particle *par = new Particle();
    par->InitializeMomentum(Particle::COORDINATE_P, Particle::COORDINATE_THETAP, p, thetap, 1.0, 1.0);
    if (eqn == "guiding-center")
        par->InitializePosition(Particle::POSITION_GUIDINGCENTER, r, 0.0, 1.0, 0.0);
    else if (eqn == "particle")
        par->InitializePosition(Particle::POSITION_PARTICLE, r, 0.0, 1.0, 0.0);

    pp->ToggleJacobianCalculation(false);

    pp->Push(par);

    delete par;
    delete conf2;
    delete conf1;
    delete mfa;

    return pp;
}

/**
 * Run the particle pusher tests.
 */
bool Test_ParticlePusher::Run(bool) {
    bool success = true;
    if (!GCMagneticAxis()) {
        this->PrintError("Exact solution to special case could not be reproduced.");
        success = false;
    } else {
        this->PrintOK("Reproduced exact solution to special case.");
    }

    if (!GuidingCenterEquation())
        success = false;

    return success;
}

