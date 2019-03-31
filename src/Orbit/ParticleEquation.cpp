/**
 * Implementation of the equations of motion
 * for a charged particle.
 */

#include <softlib/config.h>
#include <softlib/constants.h>
#include <softlib/Vector.h>
#include "Orbit/ParticleEquation.h"
#include "PhaseSpace/Particle.h"

/**
 * Classify the orbit solved by the given integrator object.
 * Since full orbits don't use 'ppar' as a coordinate, we
 * postpone the classification of the orbit (the main purpose
 * before generating the 'Orbit' object of the classification
 * is to determine whether the Jacobian was calculated properly,
 * and they are not used for full orbits).
 *
 * This way the orbit will be classified later when
 * 'ToOrbitQuantities()' is called.
 */
orbit_class_t ParticleEquation::ClassifyOrbit(Integrator<6>*) { return ORBIT_CLASS_UNKNOWN; }

/**
 * Classify the orbit based on its 'ppar' evolution.
 *
 * ppar: Parallel momentum evolution along orbit.
 * nt:   Number of time steps.
 */
orbit_class_t ParticleEquation::ClassifyOrbitPpar(const slibreal_t *ppar, const unsigned int nt) {
    unsigned int i;
    slibreal_t ppar0 = ppar[0];

    for (i = 0; i < nt; i++) {
        if (ppar[i]*ppar0 < 0)
            return ORBIT_CLASS_TRAPPED;
    }

    return ORBIT_CLASS_PASSING;
}

/**
 * Evaluate the particle velocity and acceleration
 * at time T, with particle position and velocity
 * given in zval.
 *
 * T:    Time of the values given in zval (unused).
 * zval: Particle position ([0]-[2]) and velocity ([3]-[5])
 *    in current timestep.
 * dzdt: On return, contains the calculated derivative
 *    of zval.
 * 
 * RETURNS a reference to dzdt.
 */
Vector<6>& ParticleEquation::Evaluate(const slibreal_t __UNUSED__(T), const Vector<6>& zval, Vector<6>& dzdt) {
	slibreal_t x=zval[0], y=zval[1], z=zval[2],
		px=zval[3], py=zval[4], pz=zval[5], gamma, q_gm;

	Vector<3> B = magfield->Eval(x, y, z);
	gamma = sqrt(1.0 + px*px + py*py + pz*pz);

    q_gm = particle->GetCharge() / (gamma*particle->GetMass());

	dzdt[0] = LIGHTSPEED * px / gamma;
	dzdt[1] = LIGHTSPEED * py / gamma;
	dzdt[2] = LIGHTSPEED * pz / gamma;
	dzdt[3] = q_gm * (py*B[2] - pz*B[1]);
	dzdt[4] = q_gm * (pz*B[0] - px*B[2]);
	dzdt[5] = q_gm * (px*B[1] - py*B[0]);

	return dzdt;
}

/**
 * Returns the radial particle position,
 * given a solution vector.
 *
 * x: Solution vector (6-dimensional).
 */
slibreal_t ParticleEquation::GetPositionR(
    slibreal_t x1, slibreal_t x2,
    slibreal_t __UNUSED__(x3),
    slibreal_t __UNUSED__(x4),
    slibreal_t __UNUSED__(x5),
    slibreal_t __UNUSED__(x6)
) {
    return hypot(x1, x2);
}

/**
 * Returns the vertical particle positon,
 * given a solution vector.
 *
 * x: Solution vector (6-dimensional).
 */
slibreal_t ParticleEquation::GetPositionZ(
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
 * Initialize the particle according to the
 * given specification.
 *
 * part: Particle specification.
 * zval: Vector in which to store initial value.
 *
 * RETURNS a reference to zval.
 */
Vector<6>& ParticleEquation::InitializeParticle(Particle *part, Vector<6>& zval) {
	this->particle = part;
    Vector<3> bhat, x;
    slibreal_t Babs;

    x[0] = part->GetRho();
    x[1] = 0.0;
    x[2] = part->GetZ0();

    Vector<3> B = magfield->Eval(x);
    Babs = B.Norm();
    bhat = B / Babs;

    Vector<3> p = particle->Get3Momentum(bhat);

	if (particle->GetPositionType() == Particle::POSITION_GUIDINGCENTER) {
		slibreal_t aabs, sn, cs, rhoabs;
		Vector<3> rho, pperp, a, c;

		aabs = hypot(bhat[0], bhat[1]);
		a[0] = bhat[1] / aabs;
		a[1] =-bhat[0] / aabs;
		a[2] = 0;

		c[0] = a[1]*bhat[2];
		c[1] =-a[0]*bhat[2];
		c[2] = a[0]*bhat[1] - a[1]*bhat[0];

		sn = sin(particle->GetZeta()), cs = cos(particle->GetZeta());
		rhoabs = fabs(particle->GetMass()*LIGHTSPEED*particle->GetPperp() / (particle->GetCharge() * Babs));
		rho = rhoabs * (a*cs + c*sn);
		pperp = particle->GetPperp() * (c*cs - a*sn);

		x += rho;
		p = particle->GetPpar() * bhat - pperp;
	}

	zval[0] = x[0];
	zval[1] = x[1];
	zval[2] = x[2];
	zval[3] = p[0];
	zval[4] = p[1];
	zval[5] = p[2];

	return zval;
}

/**
 * Convert a given 6D solution to this equation
 * to position (X), momentum (P), parallel momentum (ppar)
 * and perpendicular momentum (pperp).
 *
 * solution:       6D solution to this equation (1-by-(6*ntimesteps) dimensional).
 * solution2:      Secondary 6D solution to use to calculate Jacobian determinant
 *                 (set to nullptr if spatial Jacobian determinant shouldn't be calculated).
 * o:              Orbit object to store the converted result in.
 * nudge:          Nudge value used when calculating 'solution2'.
 * cl:             Orbit class (trapped, passing or unkown). If 'unknown', then this
 *                 method will try to classify the orbit.
 * forceNumerical: Force the guiding-center Jacobian to be computed numerically.
 */
void ParticleEquation::ToOrbitQuantities(
	slibreal_t *solution, slibreal_t *solution2, Orbit *o,
	slibreal_t nudge, orbit_class_t cl, bool forceNumerical
) {
    slibreal_t X,Y,Z;
    slibreal_t
        *x = o->GetX(),
        *p = o->GetP(),
        *ppar = o->GetPpar(),
        *pperp = o->GetPperp(),
        *B = o->GetB(),
        *Babs = o->GetBabs(),
        *bhat = o->GetBhat(),
        *p2 = o->GetP2(),
        *ppar2 = o->GetPpar2(),
        *pperp2 = o->GetPperp2(), 
        *gamma = o->GetGamma(),
        *gradB = o->GetGradB(),
        *curlB = o->GetCurlB(),
        **jacobianB = o->GetBJacobian();

    unsigned int i, nt = o->GetNTau();
    for (i = 0; i < nt; i++) {
        X = x[i*3+0] = solution[i*6+0];
        Y = x[i*3+1] = solution[i*6+1];
        Z = x[i*3+2] = solution[i*6+2];

        p[i*3+0] = solution[i*6+3];
        p[i*3+1] = solution[i*6+4];
        p[i*3+2] = solution[i*6+5];

        slibreal_t *_B;
        if (!o->HasBDerivatives()) {
            _B = magfield->Eval(X, Y, Z);
        } else {
            struct magnetic_field_data md = magfield->EvalDerivatives(X, Y, Z);

            _B = md.B;

            gradB[i*3+0] = md.gradB[0];
            gradB[i*3+1] = md.gradB[1];
            gradB[i*3+2] = md.gradB[2];

            curlB[i*3+0] = md.curlB[0];
            curlB[i*3+1] = md.curlB[1];
            curlB[i*3+2] = md.curlB[2];

            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    jacobianB[i*3+j][k] = md.J[j][k];
        }

        B[i*3+0] = _B[0];
        B[i*3+1] = _B[1];
        B[i*3+2] = _B[2];

        Babs[i] = sqrt(_B[0]*_B[0] + _B[1]*_B[1] + _B[2]*_B[2]);

        bhat[i*3+0] = _B[0] / Babs[i];
        bhat[i*3+1] = _B[1] / Babs[i];
        bhat[i*3+2] = _B[2] / Babs[i];

        ppar[i] = p[0]*bhat[i*3+0] + p[1]*bhat[i*3+1] + p[2]*bhat[i*3+2];

        p2[i] = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
        ppar2[i] = ppar[i]*ppar[i];
        pperp2[i] = p2[i] - ppar2[i];
        pperp[i] = sqrt(pperp2[i]);

        gamma[i] = sqrt(p2[i] + 1.0);
    }

    this->CalculateJacobians(solution, solution2, o, nudge, forceNumerical);

    if (cl == ORBIT_CLASS_UNKNOWN)
        o->SetClassification(ClassifyOrbitPpar(ppar, nt));
    else
        o->SetClassification(cl);
}

