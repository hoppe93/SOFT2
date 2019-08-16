/**
 * Module for evaluating unscreened bremsstrahlung emission in
 * the cone model.
 */

#include "Tools/Radiation/Models/Cone/ConeBremsstrahlungEmission.h"
#include "Tools/Radiation/RadiationParticle.h"

using namespace __Radiation;

ConeBremsstrahlungEmission::ConeBremsstrahlungEmission(Detector *det, MagneticField2D *mf, unsigned int nspecies, slibreal_t *Z, slibreal_t *density)
                : ConeEmission(det, mf), nspecies(nspecies)
{
    this->Z = new slibreal_t[nspecies]; 
    this->density = new slibreal_t[nspecies];
    for(unsigned int i = 0; i < nspecies; i++){
        this->Z[i] = Z[i];
        this->density[i] =  density[i];
    }
}

ConeBremsstrahlungEmission::~ConeBremsstrahlungEmission() {
    delete []this->Z; delete []this->density;
}


/**
 * Calculates the total emission and/or spectrum
 *
 * rp: Object representing particle emitting state.
 */
void ConeBremsstrahlungEmission::HandleParticle(RadiationParticle *rp, bool) {
    if (nwavelengths == 0)
        CalculateTotalEmission(rp);
    else {
        CalculateSpectrum(rp);
        IntegrateSpectrum();
    }
}

/**
 * Calculate total emission of bremsstrahlung radiation.
 * 
 * rp: Object representing the particle emitting state.
 */
void ConeBremsstrahlungEmission::CalculateTotalEmission(RadiationParticle *rp) {
    slibreal_t
        gamma = rp->GetGamma(),
        gamma2 = rp->GetGamma() * rp->GetGamma(),
        p02 = rp->GetP2(),
        p0 = sqrt(rp->GetP2()),
        logGp0 = log(gamma + p0),
        gp0 = gamma*p0;

    slibreal_t
        pf = 0,
        T1 = (12.0*gamma2 + 4.0) / (3.0*gp0) * logGp0,
        T2 =-(8.0*gamma + 6.0*p0) / (3.0*gp0*p0) * logGp0*logGp0,
        T3 =-4.0/3.0,
        T4 = 2.0 / gp0 * dilog_func(2.0*(gp0 + p02));
    for(unsigned int j = 0; j < nspecies; j++){
        pf += density[j]*Z[j]*Z[j]*this->r02Alpha;
    }

    this->power = pf * (T1 + T2 + T3 + T4);
}

/**
 * Calculate the bremsstrahlung spectrum.
 *
 * rp: Object representing the particle emitting state.
 */
void ConeBremsstrahlungEmission::CalculateSpectrum(RadiationParticle *rp) {
    slibreal_t
        gamma = rp->GetGamma(),
        gamma2 = rp->GetGamma() * rp->GetGamma(),
        p02 = rp->GetP2(),
        p0 = sqrt(rp->GetP2());
    
    slibreal_t
        pf = 0, 
        TT1 = 4.0 / 3.0,
        TT2, TT3, TT4, TT5,
        L, LT1, LT2, LT3, LT4, LT5, LTF,
        E, E2, p, p2, k, eps, eps0;
    
    for (unsigned int j = 0; j < nspecies; j++)
        pf += density[j]*Z[j]*Z[j] * this->r02Alpha;

    for (unsigned int i = 0; i < nwavelengths; i++) {
        k = wavelengths[i];

        if (k >= gamma - 1) {
            I[i] = 0;
            continue;
        }

        E = gamma - k;
        E2 = E*E;
        p2 = E2 - 1.0;
        p = sqrt(p2);

        eps  = 2.0*log(E+p);
        eps0 = 2.0*log(gamma+p0);

        TT2 =-2.0*gamma*E*(p2 + p02)/(p2*p02);
        TT3 = eps0*E/(p02*p0);
        TT4 = eps*gamma/(p2*p);
        TT5 =-eps*eps0/(p0*p);

        L = 2.0*log((gamma*E - 1.0 + p*p0)/k);

        LT1 = 8.0*gamma*E/(3.0*p0*p);
        LT2 = k*k*(gamma2*E2 + p02*p2)/(p02*p0*p2*p);

        LTF = k/(2.0*p0*p);
        LT3 = eps0 * (gamma*E + p02)/(p02*p0);
        LT4 =-eps  * (gamma*E + p2 )/(p2 *p );
        LT5 = 2.0*k*gamma*E / (p2*p02);
            
        I[i] = pf * p/(k*p0) *
            (TT1 + TT2 + TT3 + TT4 + TT5 +
             L*(LT1 + LT2 + LTF*(LT3 + LT4 + LT5)));
    }
}

/**
 * Integrates the bremsstrahlung spectrum to produce a
 * total emitted power in a given spectral range.
 */
void ConeBremsstrahlungEmission::IntegrateSpectrum() {
    unsigned int i;
    slibreal_t s = 0.5*(I[0] + I[nwavelengths-1]);

    for (i = 1; i < nwavelengths-1; i++)
        s += I[i];

    this->power = s * (wavelengths[1]-wavelengths[0]);
}

