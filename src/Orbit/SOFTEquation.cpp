/** * Implementation of basic methods in SOFTEquation.
 */

#include <omp.h>
#include <softlib/config.h>
#include <softlib/Vector.h>
#include "Orbit/SOFTEquation.h"

/**
 * Get particle radial position.
 *
 * x: Solution vector (6-dimensional).
 */
slibreal_t SOFTEquation::GetPositionR(const slibreal_t *x) {
    return GetPositionR(
        x[0], x[1], x[2],
        x[3], x[4], x[5]
    );
}
slibreal_t SOFTEquation::GetPositionR(const Vector<6>& x) {
    return GetPositionR(
        x[0], x[1], x[2],
        x[3], x[4], x[5]
    );
}

/**
 * Get particle vertical position.
 * 
 * x: Solution vector (6-dimensional).
 */
slibreal_t SOFTEquation::GetPositionZ(const slibreal_t *x) {
    return GetPositionZ(
        x[0], x[1], x[2],
        x[3], x[4], x[5]
    );
}
slibreal_t SOFTEquation::GetPositionZ(const Vector<6>& x) {
    return GetPositionZ(
        x[0], x[1], x[2],
        x[3], x[4], x[5]
    );
}

/**
 * Calculates the phase-space Jacobian denoted
 * 'J' in the SOFT documentation from the other
 * quantities of this SOFTEquation.
 *
 * solution:       Solution from first integrator.
 * solution2:      Solution from second integrator, launched
 *                 a distance 'nudge' from the first.
 * o:              Orbit object to store result in.
 * nudge:          Nudge value used when calculating 'solution2'.
 * forceNumerical: Force numerical computation of the Jacobian.
 */
void SOFTEquation::CalculateJacobians(slibreal_t *solution, slibreal_t *solution2, Orbit *o, slibreal_t nudge, bool forceNumerical) {
    unsigned int i, ti, nt = o->GetNTau();
    bool hasFlux = false;

    slibreal_t
        *p = o->GetP(),
        *p2 = o->GetP2(),
        *Babs = o->GetBabs(),
        *Beffpar = o->GetBeffpar(),
        *ppar = o->GetPpar(),
        *gamma = o->GetGamma(),
        *Jdtdrho = o->GetJdtdrho(),
        dtau = o->GetTau()[1]-o->GetTau()[0],
        J;

	// Evaluate analytical jacobian
    if (!forceNumerical && magfield->HasMagneticFlux()) {
        hasFlux = true;

        slibreal_t
            p0   = sqrt(p2[0]),
            B0   = Babs[0],
            xi0  = ppar[0] / p0,
            g0   = sqrt(1 + p0*p0),
            q    = std::abs(particle->GetCharge()),
			c    = LIGHTSPEED,
			m    = particle->GetMass(),
            *X0  = o->GetX(),
            R0   = hypot(X0[0], X0[1]),
            Z0   = X0[2],
			tauB = o->GetTau(o->GetNTau()-1),
			sigmaBp = magfield->GetCocosSigmaBp(),
			twoPi   = magfield->GetCocosEBp()==0 ? 1 : 2*M_PI;

        struct flux_diff *fd = magfield->EvalFluxDerivatives(X0);
        slibreal_t dpsi_dR    = fd->dpsi_dR;

        if (this->globset->include_drifts) {
            struct magnetic_field_data mfd = magfield->EvalDerivatives(R0, 0.0, Z0);

            slibreal_t
                bphi       = mfd.B[1] / B0,
                dB0_dR     = mfd.gradB[0],
                dbphiR0_dR = bphi + R0/B0*mfd.J[1][0] - bphi*R0/B0*dB0_dR,
				prefac = m * c*c*c*c * p0*p0*p0 * tauB / (2*M_PI*g0*q*B0);

			slibreal_t
				T1 = c*p0*bphi*R0*(1-xi0*xi0)*dB0_dR / (2*B0),
				T2 = c*p0*xi0 * dbphiR0_dR,
				T3 = (q/m)*sigmaBp/twoPi * dpsi_dR;

			J = prefac * std::abs(T1 - xi0*(T2 + T3));
        } else {
			J = c*c*c*c * p0*p0*p0 * tauB * std::abs(xi0)
				/ (2*M_PI*twoPi*g0*B0)
				* std::abs(dpsi_dR);
        }

		J *= dtau * this->particle->GetDRho();
    }

    for (i = 0; i < nt; i++) {
        // Calculate spatial Jacobian?
        if (hasFlux) {
            Jdtdrho[i] = J;
        } else if (solution2 != nullptr) {
            slibreal_t dR_dt, dZ_dt, dR_drho, dZ_drho,
                X1, Y1, Z1, X2, Y2, Z2, R1;

            ti = i*3;
            X1 = solution [i*6+0];
            Y1 = solution [i*6+1];
            Z1 = solution [i*6+2];
            X2 = solution2[i*6+0];
            Y2 = solution2[i*6+1];
            Z2 = solution2[i*6+2];
            R1 = hypot(X1, Y1);

            dR_dt = (X1*p[ti] + Y1*p[ti+1]) * LIGHTSPEED / (R1*gamma[i]) * dtau;
            dZ_dt = p[ti+2] * LIGHTSPEED / gamma[i] * dtau;

            dR_drho = hypot(X2, Y2) - R1;
            dZ_drho = Z2 - Z1;

            Jdtdrho[i] = R1 * fabs(dR_dt*dZ_drho - dZ_dt*dR_drho) / nudge;
            if (this->particle != 0)
                Jdtdrho[i] *= this->particle->GetDRho();

            // Jp
            Jdtdrho[i] *= Beffpar[i] / Babs[0] * fabs(ppar[0] / ppar[i]);
        }
    }
}

