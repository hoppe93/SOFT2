/**
 * Implementation of the relativistic guiding-center
 * equations of motion.
 *
 * Things to test:
 * - If guiding-centers are correctly placed when the
 *   particle position is given.
 */

#include <softlib/config.h>
#include <softlib/constants.h>
#include "Orbit/GuidingCenterEquation.h"
#include "Orbit/ParticlePusher.h"

/**
 * Classify the orbit solved by the given integrator object.
 *
 * intg: Integrator object containing the orbit solution.
 */
orbit_class_t GuidingCenterEquation::ClassifyOrbit(Integrator<6> *intg) {
    unsigned int nsteps = intg->StepsTaken(), i;
    slibreal_t *sol = intg->SolutionAt(0), ppar0;

    ppar0 = sol[COORD_PPAR];
    for (i = 0; i < nsteps; i++) {
        if (sol[i*6 + COORD_PPAR]*ppar0 < 0)
            return ORBIT_CLASS_TRAPPED;
    }

    return ORBIT_CLASS_PASSING;
}

/**
 * Classify the orbit based on its 'ppar' evolution.
 *
 * ppar: Parallel momentum evolution along orbit.
 * nt:   Number of time steps.
 */
orbit_class_t GuidingCenterEquation::ClassifyOrbitPpar(const slibreal_t *ppar, const unsigned int nt) {
    unsigned int i;
    slibreal_t ppar0 = ppar[0];

    for (i = 0; i < nt; i++) {
        if (ppar[i]*ppar0 < 0)
            return ORBIT_CLASS_TRAPPED;
    }

    return ORBIT_CLASS_PASSING;
}

/**
 * Evaluate the guiding-center equations of motion.
 * 
 * T:    Current time.
 * zval: Solution at time T.
 * dzdt: Contains derivative of zval with respect to time
 *       at time T on return.
 *
 * RETURNS dzdt, which has components
 *   [0]: Guiding-center velocity in x-direction (X-dot).
 *   [1]: Guiding-center velocity in y-direction (Y-dot).
 *   [2]: Guiding-center velocity in z-direction (Z-dot).
 *   [3]: Guiding-center parallel force (ppar-dot).
 *   [4]: Time-rate of change of particle magnetic moment (mu-dot = 0).
 *   [5]: Particle gyration frequency (zeta-dot).
 */
Vector<6>& GuidingCenterEquation::Evaluate(const slibreal_t __UNUSED__(T), const Vector<6>& zval, Vector<6>& dzdt, slibreal_t *gamma) {
	slibreal_t
        x = zval[COORD_X],
        y = zval[COORD_Y],
        z = zval[COORD_Z],
		ppar = zval[COORD_PPAR],
        mu = zval[COORD_MU];
    slibreal_t
        c = LIGHTSPEED,
        m = this->particle->GetMass(),
        q = this->particle->GetCharge(),
        Bstarpar, p2, gmm;

    struct magnetic_field_data mfd = magfield->EvalDerivatives(x, y, z);
    Vector<3> Bstar, B(mfd.B), curlB(mfd.curlB), curlBhat, gradB(mfd.gradB);
    Vector<3> bhat(B/mfd.Babs), Xdot, bhatXgradB;

    p2 = ppar*ppar + 2.0*mfd.Babs*mu/(m*c*c);
    gmm = sqrt(p2+1.0);

    // Calculate B* and Xdot
    if (this->include_drifts) {
        Vector<3>::Cross(bhat, gradB, bhatXgradB);
        curlBhat = (curlB + bhatXgradB) / mfd.Babs;
        Bstar = B + (m*c*ppar / q) * curlBhat;
        Bstarpar = bhat.Dot(Bstar);

        Xdot = (c*ppar*Bstar + (mu/q)*bhatXgradB) / (gmm*Bstarpar);
    } else {
        Bstar = B;
        Bstarpar = mfd.Babs;

        Xdot = (c*ppar/(gmm*Bstarpar))*Bstar;
    }

    dzdt[COORD_X]    = Xdot[0];
    dzdt[COORD_Y]    = Xdot[1];
    dzdt[COORD_Z]    = Xdot[2];
    dzdt[COORD_PPAR] =-mu/(gmm*m*c*Bstarpar) * gradB.Dot(Bstar);
    dzdt[COORD_MU]   = 0.0;
    dzdt[COORD_ZETA] = q*mfd.Babs / (gmm*m);

    if (gamma != nullptr)
        *gamma = gmm;

    return dzdt;
}

/**
 * Returns the radial coordinate of the guiding-center,
 * as defined by the given 6-dimensional solution vector.
 *
 * x: X coordinate of guiding-center.
 * y: Y coordinate of guiding-center.
 *
 * All other parameters are unused.
 */
slibreal_t GuidingCenterEquation::GetPositionR(
    slibreal_t x, slibreal_t y,
    slibreal_t __UNUSED__(x3),
    slibreal_t __UNUSED__(x4),
    slibreal_t __UNUSED__(x5),
    slibreal_t __UNUSED__(x6)
) {
    return hypot(x, y);
}

/**
 * Returns the vertical coordinate of the guiding-center,
 * as defined by the given 6-dimensional solution vector.
 *
 * x1: Unused.
 * x2: Unused.
 * z:  Z (vertical) coordinate of the guiding-center.
 *
 * All other parameters are unused.
 */
slibreal_t GuidingCenterEquation::GetPositionZ(
    slibreal_t __UNUSED__(x1),
    slibreal_t __UNUSED__(x2),
    slibreal_t z,
    slibreal_t __UNUSED__(x4),
    slibreal_t __UNUSED__(x5),
    slibreal_t __UNUSED__(x6)
) {
    return z;
}

/**
 * Initialize the guiding-center equation with
 * the given particle object.
 *
 * part: Particle specifying settings to use.
 * zval: Solution object to initialize.
 *
 * RETURNS a reference to zval.
 */
Vector<6>& GuidingCenterEquation::InitializeParticle(Particle *part, Vector<6>& zval) {
    this->particle = part;
    Vector<3> x;
    slibreal_t Babs, mu, pp;

    x[0] = part->GetRho();
    x[1] = 0.0;
    x[2] = part->GetZ0();

    Vector<3> B = magfield->Eval(x);
    Babs = B.Norm();
    Vector<3> bhat = B / Babs;

    slibreal_t
        m = part->GetMass(),
        c = LIGHTSPEED;

    pp = part->GetPperp();
    mu = pp*pp*m*c*c / (2.0*Babs);

    if (particle->GetPositionType() == Particle::POSITION_PARTICLE) {
        Vector<3> rho, pperp;

        pperp = part->Get3Momentum(bhat) - part->GetPpar() * bhat;

        rho[0] = bhat[2]*pperp[1] - bhat[1]*pperp[2];
        rho[1] = bhat[0]*pperp[2] - bhat[2]*pperp[0];
        rho[2] = bhat[1]*pperp[0] - bhat[0]*pperp[1];

        rho /= part->GetCharge()*Babs / (m*c);

        x += rho;
    }

    zval[COORD_X]    = x[0];
    zval[COORD_Y]    = x[1];
    zval[COORD_Z]    = x[2];
    zval[COORD_PPAR] = part->GetPpar();
    zval[COORD_MU]   = mu;
    zval[COORD_ZETA] = part->GetZeta();

    return zval;
}

/**
 * Convert a given 6D solution to this equation
 * to an Orbit object.
 *
 * solution:       6D solution to this equation (1-by-(6*ntimesteps) dimensional).
 * solution2:      Secondary 6D solution to use to calculate Jacobian determinant
 *                 (set to nullptr if spatial Jacobian determinant shouldn't be calculated).
 * timingSolution: Solution used for finding the poloidal transit time when solving a
 *                 particle equation motion. Hence not used here.
 * o:              Orbit object to store converted result in.
 * nudge:          Nudge value used to calculate 'solution2'.
 * cl:             Orbit classification (trapped, passing or unknown). If 'unknown', this
 *                 method will try to classify the orbit.
 * forceNumerical: Force the guiding-center Jacobian to be computed numerically.
 */
void GuidingCenterEquation::ToOrbitQuantities(
	slibreal_t *solution, slibreal_t *solution2, slibreal_t*,
    Orbit *o, slibreal_t nudge, orbit_class_t cl, bool forceNumerical
) {
    Vector<6> dzdt;
    slibreal_t X,Y,Z;
    slibreal_t
        *x = o->GetX(),
        *p = o->GetP(),
        *ppar = o->GetPpar(),
        *pperp = o->GetPperp(),
        *B = o->GetB(),
        *Babs = o->GetBabs(),
        *Beffpar = o->GetBeffpar(),
        *bhat = o->GetBhat(),
        *p2 = o->GetP2(),
        *ppar2 = o->GetPpar2(),
        *pperp2 = o->GetPperp2(), 
        *gamma = o->GetGamma(),
        *gradB = o->GetGradB(),
        *curlB = o->GetCurlB(),
        **jacobianB = o->GetBJacobian();

    slibreal_t
        c = LIGHTSPEED,
        m = this->particle->GetMass(),
        q = this->particle->GetCharge();

    unsigned int i, nt = o->GetNTau();
    for (i = 0; i < nt; i++) {
        X = x[i*3+0] = solution[i*6+COORD_X];
        Y = x[i*3+1] = solution[i*6+COORD_Y];
        Z = x[i*3+2] = solution[i*6+COORD_Z];

        Vector<6> z(solution+(i*6));
        Evaluate(0.0, z, dzdt, gamma+i);
        p[i*3+0] = dzdt[0] * gamma[i] / LIGHTSPEED;
        p[i*3+1] = dzdt[1] * gamma[i] / LIGHTSPEED;
        p[i*3+2] = dzdt[2] * gamma[i] / LIGHTSPEED;

        ppar[i] = solution[i*6+COORD_PPAR];

        if (!o->HasBDerivatives()) {
            slibreal_t *_B = magfield->Eval(X, Y, Z);

            B[i*3+0] = _B[0];
            B[i*3+1] = _B[1];
            B[i*3+2] = _B[2];

            Babs[i] = sqrt(_B[0]*_B[0] + _B[1]*_B[1] + _B[2]*_B[2]);

            Beffpar[i] = Babs[i];
        } else {
            struct magnetic_field_data md = magfield->EvalDerivatives(X, Y, Z);
            slibreal_t *_B = md.B;

            gradB[i*3+0] = md.gradB[0];
            gradB[i*3+1] = md.gradB[1];
            gradB[i*3+2] = md.gradB[2];

            curlB[i*3+0] = md.curlB[0];
            curlB[i*3+1] = md.curlB[1];
            curlB[i*3+2] = md.curlB[2];

            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    jacobianB[i*3+j][k] = md.J[j][k];

            B[i*3+0] = _B[0];
            B[i*3+1] = _B[1];
            B[i*3+2] = _B[2];

            Babs[i] = sqrt(_B[0]*_B[0] + _B[1]*_B[1] + _B[2]*_B[2]);

            bhat[i*3+0] = _B[0] / Babs[i];
            bhat[i*3+1] = _B[1] / Babs[i];
            bhat[i*3+2] = _B[2] / Babs[i];

            Vector<3> _bhat(bhat+3*i), _gradB(gradB+3*i), _curlB(curlB+3*i), __B(_B);

            Vector<3> bhatXgradB = Vector<3>::Cross(_bhat, _gradB);
            Vector<3> curlBhat = (_curlB + bhatXgradB) / Babs[i];
            Vector<3> Bstar = __B + (m*c*ppar[i] / q) * curlBhat;

            Beffpar[i] = _bhat.Dot(Bstar);
        }

        ppar2[i] = ppar[i]*ppar[i];
        pperp2[i] = solution[i*6+COORD_MU] * 2.0*Babs[i] / (o->GetMass()*LIGHTSPEED*LIGHTSPEED);
        p2[i] = ppar2[i] + pperp2[i];
        pperp[i] = sqrt(pperp2[i]);
    }

    this->CalculateJacobians(solution, solution2, o, nudge, forceNumerical);

    if (cl == ORBIT_CLASS_UNKNOWN)
        o->SetClassification(ClassifyOrbitPpar(ppar, nt));
    else
        o->SetClassification(cl);
}

/**
 * Enable or disable the guiding-center drift
 * velocity terms.
 *
 * enable: If true, adds the guiding-center drift velocity
 *   terms to the equations of motion. Otherwise only
 *   guiding-center motino along magnetic field lines will
 *   be considered.
 */
void GuidingCenterEquation::ToggleDrifts(bool enable) { this->include_drifts = enable; }

