/**
 * Implementation of momentum in FOGC.
 */

#include <softlib/config.h>
#include <softlib/constants.h>
#include <softlib/Vector.h>

#include "GOAMomentum.h"

/**
 * Constructor.
 *
 * X:  Guiding-center position.
 * mf: Magnetic field to work in.
 * q:  Particle charge (in units of the elementary charge).
 */
GOAMomentum::GOAMomentum(Vector<3>& X, MagneticField2D *mf, slibreal_t q) {
    SetGCPosition(X, mf);
    this->charge = q * ELECTRON_CHARGE;
}

/**
 * Evaluate the particle momentum up to and including
 * the first order corrections in guiding-center theory.
 *
 * X:     Guiding-center position.
 * ppar:  Particle parallel momentum.
 * pperp: Particle perpendicular momentum.
 * zeta:  Particle gyrophase.
 * p:     Vector containing the particle momentum on return.
 */
Vector<3>& GOAMomentum::ParticleMomentum(
    slibreal_t ppar, slibreal_t pperp, slibreal_t zeta,
    Vector<3>& P, Vector<3>& p
) {
    Vector<3> p1;

    // Guiding-center momentum P
    P = GuidingCenterMomentum(ppar, pperp, P);

    p = P + Momentum1(ppar, pperp, zeta, p1);

    return p;
}

/**
 * Evaluate the guiding-center drift momentum,
 * p_{drift}.
 *
 * ppar:     Particle parallel momentum.
 * pperp:    Particle perpendicular momentum.
 * pdrift:   Contains guiding-center drift-momentum on return.
 */
Vector<3>& GOAMomentum::GuidingCenterMomentum(
    slibreal_t ppar, slibreal_t pperp, Vector<3>& P
) {
    Vector<3> Beff, muterm, parterm;
    slibreal_t Beffpar, q = this->charge;

    Beff = B + ppar*ELECTRON_MASS*LIGHTSPEED/q*curlBhat;
    Beffpar = bhat.Dot(Beff);

    Vector<3>::Cross(gradB, bhat, muterm);
    muterm = muterm * (ELECTRON_MASS*LIGHTSPEED*pperp*pperp/(2.0*q*Beffpar*Babs));
    parterm = (ppar/Beffpar)*Beff;
    P = parterm - muterm;

    return P;
}

/**
 * Constructs a coordinate system (1, 2, b), such that
 * the unit vectors 1, 2 & b are mutually orthogonal.
 */
void GOAMomentum::InitLocalCoordinateSystem() {
    if (fabs(bhat[2]) > fabs(bhat[1]) && fabs(bhat[2]) > fabs(bhat[0])) {
        one[0] = 1.0;
        one[1] = 1.0;
        one[2] = -(bhat[0] + bhat[1]) / bhat[2];
    } else if (fabs(bhat[1]) > fabs(bhat[0])) {
        one[0] = 1.0;
        one[1] = -(bhat[2] + bhat[0]) / bhat[1];
        one[2] = 1.0;
    } else {
        one[0] = -(bhat[1] + bhat[2]) / bhat[0];
        one[1] = 1.0;
        one[2] = 1.0;
    }

    one.Normalize();
    Vector<3>::Cross(bhat, one, two);
}

/**
 * First order correction to the particle momentum.
 * Note that also the zeroth-order term, pperp*perphat,
 * is included here.
 *
 * ppar:  Particle parallel momentum.
 * pperp: Particle perpendicular momentum.
 * zeta:  Particle gyrophase.
 * p1:    Vector to store output in.
 */
Vector<3>& GOAMomentum::Momentum1(
    slibreal_t ppar, slibreal_t pperp,
    slibreal_t zeta, Vector<3>& p1
) {
    slibreal_t c = LIGHTSPEED,
               m = ELECTRON_MASS,
               q = this->charge;
    slibreal_t rL = m*c*pperp / (fabs(q)*Babs),
               //mu = m*c*c*pperp*pperp / (2.0*Babs);
               mu = m*c*c*pperp*pperp / (2.0*Babs);

    Vector<3> rhohat  = one*cos(zeta) - two*sin(zeta),
              perphat =-one*sin(zeta) - two*cos(zeta);
    slibreal_t bterm, pterm, rterm, a1GradB, a2GradB;

    Vector<3> gradLogB(gradB / Babs);

    GetContractions(rhohat, perphat, &a1GradB, &a2GradB);

    // Collect terms
    slibreal_t bt1, bt2,
               pt1, pt2, pt3,
               rt1, rt2, rt3;
    /*bterm = (2.0*m*c*ppar*pperp/(q*Babs) - rL*ppar)*kappa.Dot(rhohat) + 0.5*rL*pperp*a1GradB;
    pterm = q*Babs*rL/(m*c) + 0.5*rL*pperp*rhohat.Dot(gradLogB) - 0.5*m*c*ppar*pperp/(q*Babs)*a1GradB;
    rterm = -m*c*ppar*pperp/(q*Babs)*a2GradB + 0.5*rL*pperp*perphat.Dot(gradLogB) - 0.5*ppar*rL*bhat.Dot(gradLogB);*/

    // OLD IMPLEMENTATION BASED ON OLA'S NOTES
    /*
    bt1 = (2.0*m*c*ppar*pperp/(q*Babs) - rL*ppar)*kappa.Dot(rhohat);
    bt2 = 0.5*rL*pperp*a1GradB;

    pt1 = q*Babs*rL/(m*c);
    pt2 = 0.5*rL*pperp*rhohat.Dot(gradLogB);
    pt3 =-0.5*m*c*ppar*pperp/(q*Babs)*a1GradB;

    rt1 = -m*c*ppar*pperp/(q*Babs)*a2GradB;
    rt2 = 0.5*rL*pperp*perphat.Dot(gradLogB);
    rt3 =-0.5*ppar*rL*bhat.Dot(gradLogB);
    */

    // NEW IMPLEMENTATION, BASED ON MY DERIVATION
    bt1 = ppar*rL*kappa.Dot(rhohat);
    bt2 = mu * a1GradB / (c*q);

    pt1 = q*Babs*rL/(m*c);
    pt2 = 0.5*q*Babs*rL*rL * rhohat.Dot(gradLogB) / (m*c);
    pt3 =-0.5*rL*ppar*a1GradB;

    rt1 = 0.5*q*Babs*rL*rL*perphat.Dot(gradLogB) / (m*c);
    rt2 =-0.5*ppar*rL*bhat.Dot(gradLogB);
    rt3 =-rL*ppar*a2GradB;

    bterm = bt1 + bt2;
    pterm = pt1 + pt2 + pt3;
    rterm = rt1 + rt2 + rt3;

    p1 = bterm*bhat + pterm*perphat + rterm*rhohat;
    return p1;
}

/**
 * Calculate the gyro radius vector rho.
 *
 * ppar:  Particle parallel momentum.
 * pperp: Particle perpendicular momentum.
 * zeta:  Particle gyrophase.
 * rho:   Vector to store output in.
 */
Vector<3>& GOAMomentum::GyroVector(
    slibreal_t ppar, slibreal_t pperp,
    slibreal_t zeta, Vector<3>& rho
) {
    slibreal_t rL = ELECTRON_MASS*LIGHTSPEED*pperp / (ELECTRON_CHARGE*Babs),
        rhopar = ELECTRON_MASS*LIGHTSPEED*ppar / (ELECTRON_CHARGE*Babs);

    Vector<3> rhohat  = one*cos(zeta) - two*sin(zeta),
              perphat =-one*sin(zeta) - two*cos(zeta);
    Vector<3> rho0, rho1, dRho0Dzeta, dRho0Dmu;
    Vector<3> gradLogB(gradB / Babs);

    slibreal_t bterm, rhoterm, perpterm;
    slibreal_t a1GradB, a2GradB, tauB;

    GetContractions(rhohat, perphat, &a1GradB, &a2GradB);

    dRho0Dzeta = rL*perphat;
    dRho0Dmu   = rhohat/(LIGHTSPEED*ELECTRON_CHARGE*pperp);
    tauB       = bhat.Dot(curlBhat);

    // Zeroth-order
    rho0 = rL * rhohat;

    // First-order
    bterm    = 2.0*rhopar*dRho0Dzeta.Dot(kappa) + 0.5*rL*rL*(a2GradB + 0.5*bhat.Dot(gradB)/Babs);
    rhoterm  = 0.5*rhopar*(tauB - a1GradB) - 0.5*rL*rhohat.Dot(gradLogB + (2.0*ppar*ppar/(pperp*pperp))*kappa);
    perpterm = rhopar*a2GradB - dRho0Dzeta.Dot(gradLogB + (ppar*ppar/(pperp*pperp))*kappa);

    rho1 = bterm*bhat + rhoterm*rL*rhohat + perpterm*rL*perphat;
    rho = rho0 + rho1;

    return rho;
}

/* Direct implementation of Eero's equations.
 * Opposite sign to the above.
Vector<3>& Momentum::GyroVector(
    slibreal_t ppar, slibreal_t pperp,
    slibreal_t zeta, Vector<3>& rho
) {
    slibreal_t rL = ELECTRON_MASS*LIGHTSPEED*pperp / (ELECTRON_CHARGE*Babs),
        rhopar = ELECTRON_MASS*LIGHTSPEED*ppar / (ELECTRON_CHARGE*Babs),
        mu = ELECTRON_MASS*LIGHTSPEED*LIGHTSPEED*pperp*pperp / (2.0*Babs);

    Vector<3> rhohat  = one*cos(zeta) - two*sin(zeta),
              perphat =-one*sin(zeta) - two*cos(zeta);
    Vector<3> rho0, rho1, dRho0Dzeta, dRho0Dmu;
    Vector<3> gradLogB(gradB / Babs);

    slibreal_t bterm, rhoterm, perpterm;
    slibreal_t a1GradB, a2GradB, tauB, t1, t2, t3;

    GetContractions(rhohat, perphat, &a1GradB, &a2GradB);

    dRho0Dzeta = rL*perphat;
    dRho0Dmu   = rhohat/(LIGHTSPEED*ELECTRON_CHARGE*pperp);
    tauB       = bhat.Dot(curlBhat);

    // Zeroth-order
    rho0 = rL * rhohat;

    // First-order
    bterm    = 2.0*rhopar*kappa.Dot(dRho0Dzeta) - ELECTRON_MASS*mu/(ELECTRON_CHARGE*ELECTRON_CHARGE*Babs)*(a2GradB + 0.5*bhat.Dot(gradB)/Babs);
    rhoterm  = 0.5*rhopar*(tauB - a1GradB) + dRho0Dmu.Dot(mu*gradLogB + (ELECTRON_MASS*LIGHTSPEED*LIGHTSPEED*ppar*ppar/Babs) * kappa);
    perpterm = rhopar*a2GradB + dRho0Dzeta.Dot(gradLogB + (ELECTRON_MASS*LIGHTSPEED*LIGHTSPEED*ppar*ppar / (2.0*mu*Babs)) * kappa);

    rho1 = -bterm*bhat - rhoterm*rho0 - perpterm*dRho0Dzeta;
    rho = rho0 + rho1;

    return rho;
}*/

void GOAMomentum::GetContractions(
    Vector<3>& rhohat, Vector<3>& perphat,
    slibreal_t *a1GradB, slibreal_t *a2GradB
) {
    // Gradients of individual components of B
    Vector<3> dBx(J[0]), dBy(J[1]), dBz(J[2]);
    // Corresponding gradients of individual components of bhat
    Vector<3> dbx = (dBx - bhat[0]*gradB) / Babs,
              dby = (dBy - bhat[1]*gradB) / Babs,
              dbz = (dBz - bhat[2]*gradB) / Babs;

    // Tensor products
                    // rhohat . (perphat . nabla) bhat
    *a1GradB=-0.5*(rhohat[0]*perphat.Dot(dbx) +
                   rhohat[1]*perphat.Dot(dby) +
                   rhohat[2]*perphat.Dot(dbz) +
                   // perphat . (rhohat . nabla) bhat
                   perphat[0]*rhohat.Dot(dbx) +
                   perphat[1]*rhohat.Dot(dby) +
                   perphat[2]*rhohat.Dot(dbz));

                    // perphat . (perphat . nabla) bhat
    *a2GradB= 0.25*(perphat[0]*perphat.Dot(dbx) +
                   perphat[1]*perphat.Dot(dby) +
                   perphat[2]*perphat.Dot(dbz) -
                   // rhohat . (rhohat . nabla) bhat
                   rhohat[0]*rhohat.Dot(dbx) -
                   rhohat[1]*rhohat.Dot(dby) -
                   rhohat[2]*rhohat.Dot(dbz));
}

/**
 * Calculates the particle acceleration induced by
 * the magnetic field to first order (i.e. beta-dot).
 *
 * p:       Particle momentum.
 * ppar:    Particle parallel momentum.
 * pperp:   Particle perpendicular momentum.
 * zeta:    Gyrophase.
 * betadot: Contains acceleration on return.
 */
Vector<3>& GOAMomentum::Acceleration(
    Vector<3>& p, slibreal_t ppar, slibreal_t pperp, slibreal_t zeta,
    Vector<3>& betadot
) {
    slibreal_t rL = ELECTRON_MASS*LIGHTSPEED*pperp / (fabs(this->charge)*Babs),
        prefac, gamma = sqrt(p.Dot(p) + 1.0);
    Vector<3> t1, t2, t3, gradBrho0;
    Vector<3> rhohat  = one*cos(zeta) - two*sin(zeta),
              perphat =-one*sin(zeta) - two*cos(zeta);
    Vector<3> rho0 = rL * rhohat;
    Vector<3> dBx(J[0]), dBy(J[1]), dBz(J[2]);
    Vector<3> part, perpt;

    gradBrho0[0] = rho0.Dot(dBx);
    gradBrho0[1] = rho0.Dot(dBy);
    gradBrho0[2] = rho0.Dot(dBz);

    part = ppar*bhat;
    perpt = pperp*perphat;

    Vector<3>::Cross(p, B, t1);
    Vector<3>::Cross(part, gradBrho0, t2);
    Vector<3>::Cross(perpt, gradBrho0, t3);

    prefac = this->charge/(ELECTRON_MASS*gamma);

    betadot = prefac * (t1 + t2 + t3);
    return betadot;
}

/**
 * Set the guiding-center position and calculate
 * related magnetic field quantities.
 *
 * X:  New guiding-center position.
 * mf: Magnetic field.
 */
void GOAMomentum::SetGCPosition(Vector<3>& X, MagneticField2D *magfield) {
    struct magnetic_field_data mfd = magfield->EvalDerivatives(X);
    Vector<3> bGradB;

    gradB = Vector<3>(mfd.gradB);
    B = Vector<3>(mfd.B);
    curlB = Vector<3>(mfd.curlB);
    Babs = mfd.Babs;

    bhat = B / mfd.Babs;
    curlBhat = (curlB + Vector<3>::Cross(bhat, gradB, bGradB)) / Babs;

    // kappa = (bhat . nabla) bhat =
    //       = -bhat x (curl bhat) =
    //       = (curl bhat) x bhat.
    Vector<3>::Cross(curlBhat, bhat, kappa);
    
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            J[i][j] = mfd.J[i][j];
        }
    }

    InitLocalCoordinateSystem();
}

/**
 * Print all magnetic field quantities to stdout.
 */
void GOAMomentum::DumpMF() {
    printf("|B|      = %12.8e\n", Babs);
    printf("B        = %12.8e, %12.8e, %12.8e\n", B[0], B[1], B[2]);
    printf("bhat     = %12.8e, %12.8e, %12.8e\n", bhat[0], bhat[1], bhat[2]);
    printf("gradB    = %12.8e, %12.8e, %12.8e\n", gradB[0], gradB[1], gradB[2]);
    printf("curlB    = %12.8e, %12.8e, %12.8e\n", curlB[0], curlB[1], curlB[2]);
    printf("curlBhat = %12.8e, %12.8e, %12.8e\n", curlBhat[0], curlBhat[1], curlBhat[2]);
    printf("kappa    = %12.8e, %12.8e, %12.8e\n", kappa[0], kappa[1], kappa[2]);
    printf("J        = %12.8e, %12.8e, %12.8e\n"
           "           %12.8e, %12.8e, %12.8e\n"
           "           %12.8e, %12.8e, %12.8e\n",
           J[0][0], J[0][1], J[0][2],
           J[1][0], J[1][1], J[1][2],
           J[2][0], J[2][1], J[2][2]);
}

