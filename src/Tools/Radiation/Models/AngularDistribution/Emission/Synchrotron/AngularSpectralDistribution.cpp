/**
 * Implementation of the angular and spectral distribution
 * of synchrotron radiation, including polarization components.
 */

#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"
#include "Tools/Radiation/synchrotron_func.h"

using namespace std;
using namespace __Radiation;

/*************************
 * PREPARATION FUNCTIONS *
 *************************/
void ADSynchrotronEmission::PreparePolarization(RadiationParticle *rp) {
    PrepareSpectrum(rp);

    this->Efield.nE = this->nwavelengths;

    this->Efield.Ex2  = new slibreal_t[this->Efield.nE];
    this->Efield.Ey2  = new slibreal_t[this->Efield.nE];
    this->Efield.ExEy = new slibreal_t[this->Efield.nE];
}
void ADSynchrotronEmission::PrepareSpectrum(RadiationParticle *rp) {
    this->gamma = rp->GetGamma();
    this->igamma2 = 1.0 / (gamma*gamma);
    this->gamma3 = gamma*gamma*gamma;

    slibreal_t
        m = rp->GetMass(),
        q = fabs(rp->GetCharge()),
        B = rp->GetB(),
        p2 = rp->GetP2(),
        p = sqrt(p2),
        beta2 = p2 / (1.0 + p2),
        betapar2 = rp->GetPpar()*rp->GetPpar() / (1.0 + p2),
        gammapar = 1.0 / sqrt(1.0 - betapar2);

    this->prefactor = q*q*q*beta2*B / (16.0*M_PI*EPS0*gamma*m);
    this->beta = sqrt(beta2);
    this->cosThetap = rp->GetPpar() / p;
    this->sinThetap = rp->GetPperp() / p;
    this->betapar = beta*cosThetap;
	this->Vsign = (rp->GetP().Dot(rp->GetBvec())) >= 0 ? 1 : -1;

    this->lambdac = 4.0*M_PI*LIGHTSPEED*gammapar*igamma2*m / (3.0*q*B);
}

/**
 * Wrapper for calculating polarization of
 * synchrotron radiation.
 */
void ADSynchrotronEmission::CalculatePolarization(
    RadiationParticle *rp,
    Vector<3> &n, slibreal_t sinMu,  slibreal_t cosMu
) { __CalculateSpectrum<true>(rp, n, sinMu, cosMu); }

/**
 * Wrapper for calculating angular and spectral
 * distribution of synchrotron radiation.
 */
void ADSynchrotronEmission::CalculateSpectrum(
    Vector<3> &n, slibreal_t sinMu,  slibreal_t cosMu
) { __CalculateSpectrum<false>(nullptr, n, sinMu, cosMu); }

/**
 * Function for calculating the spectrum of emitted synchrotron
 * radiation, including its polarization components.
 * 
 * calculatePolarization: If true, calculates the polarization
 *                        components of synchrotron radiation.
 * rp:                    Object representing the particle
 *                        emitting state.
 */
template<bool calculatePolarization>
void ADSynchrotronEmission::__CalculateSpectrum(
    RadiationParticle *rp,
    Vector<3>&, slibreal_t sinMu,  slibreal_t cosMu
) {
    unsigned int i;
    slibreal_t
        cosPsi = (Vsign*cosMu*cosThetap + sinMu*sinThetap),
        sinPsi2 = 1.0 - cosPsi*cosPsi,
        mcospsi = 1.0-beta*cosPsi,
        xifac = gamma3*lambdac*sqrt(mcospsi*mcospsi*mcospsi/(0.5*beta*cosPsi));

    if (cosPsi < 0) {
        for (i = 0; i < nwavelengths; i++)
            I[i] = 0.0;
        return;
    }

    if (calculatePolarization) {
        Vector<3> bhat = this->magfield->Eval(rp->GetPosition());
        bhat.Normalize();

        slibreal_t r = rp->GetRCP().Norm();
        Vector<3> rcphat = rp->GetRCP() / r;
        slibreal_t nDotB = rcphat.Dot(bhat);

        this->Efield.zhat = rcphat;
        this->Efield.yhat =-Vector<3>::Cross(rcphat, bhat) / (slibreal_t)sqrt(1.0 - nDotB*nDotB);
        this->Efield.xhat = Vector<3>::Cross(this->Efield.yhat, this->Efield.zhat);
    }

    for (i = 0; i < nwavelengths; i++) {
        slibreal_t
            l = wavelengths[i],
            pfac = this->prefactor/(l*l * beta*cosPsi*(1.0 - beta*cosPsi)) * (1 - betapar*cosMu),
            //pfac = 1.0/(l*l * beta*cosPsi*(1.0 - beta*cosPsi)) * (1 - betapar*cosMu),
            fac13 = (0.5*beta*cosPsi*sinPsi2)/(1.0 - beta*cosPsi),
            xi = xifac / l;

        slibreal_t
            xK13 = synchrotron_func3(xi),
            xK23 = synchrotron_func2(xi);

        if (calculatePolarization) {
            Efield.Ex2[i] = pfac*xK23*xK23;
            Efield.Ey2[i] = pfac*fac13*xK13*xK13;
            Efield.ExEy[i] = sqrt(fac13)*pfac*xK13*xK23;
        } else
            I[i] = pfac * (xK23*xK23 + fac13*xK13*xK13);
    }

    if (calculatePolarization)
        this->detector->GetOptics()->ApplyOptics(Efield, I, Q, U, V);
}

