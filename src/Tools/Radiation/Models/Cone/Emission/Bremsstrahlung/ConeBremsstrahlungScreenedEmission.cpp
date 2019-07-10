/**
 * Module for evaluating screened bremsstrahlung emission in
 * the cone model.
 */

#include "Tools/Radiation/Models/Cone/ConeBremsstrahlungScreenedEmission.h"
#include "Tools/Radiation/RadiationParticle.h"

using namespace __Radiation;
/* ConeBremsstrahlungScreenedEmission(Detector *det, MagneticField2D *mf, unsigned int nspecies, slibreal_t *Z, slibreal_t *Z0, slibreal_t *density)
                : ConeEmission(det, mf), nspecies(nspecies), Z(Z), Z0(Z0), density(density) 
{
    //this->r02 = ELECTRON_CHARGE*ELECTRON_CHARGE*ELECTRON_CHARGE*ELECTRON_CHARGE/(16*M_PI*M_PI*EPS0*EPS0*LIGHTSPEED*LIGHTSPEED*LIGHTSPEED*LIGHTSPEED*ELECTRON_MASS*ELECTRON_MASS);
            }

/* //This part does not currently compile, might be source of problem

/**
 * Calculates the total emission and/or spectrum and/or
 * Stokes parameters of synchrotron radiation.
 *
 * rp:           Object representing particle emitting state.
 * 
 */
void ConeBremsstrahlungScreenedEmission::HandleParticle(RadiationParticle *rp, bool) {
    if (nwavelengths == 0)
        CalculateTotalEmission();
    else {
        CalculateSpectrum(rp);
        IntegrateSpectrum();
    }
}

/**
 * Calculate total emission of bremsstrahlung, heavily approximated.
 * Uses eg. 4BS from Koch and Motz
 * 
 *
 */
void ConeBremsstrahlungScreenedEmission::CalculateTotalEmission() { //Test will not work with RadiationParticle as argument here
    slibreal_t I = 0;   
    //this->r02 = ELECTRON_CHARGE*ELECTRON_CHARGE*ELECTRON_CHARGE*ELECTRON_CHARGE/(16*M_PI*M_PI*EPS0*EPS0*LIGHTSPEED*LIGHTSPEED*LIGHTSPEED*LIGHTSPEED*ELECTRON_MASS*ELECTRON_MASS);

    for(unsigned int i=0; i < nspecies; i++){
        I = I + density[i]*Calculate4BS(Z[i]);
    }
    this->power = I;
}

/**
 * Calculate eq. 4BS (Koch & Motz) for a given species
 */

slibreal_t ConeBremsstrahlungScreenedEmission::Calculate4BS(slibreal_t Z) {
    slibreal_t Zfakt = 4*Z*Z*r02/137;
    slibreal_t lnfakt = log(183) + 0.5*log(Z) + 0.0555555555555556;
    
    return Zfakt*lnfakt;
}

/**
 * Calculate the bremsstrahlung spectrum. NOTE!! HAS NOT HAD BEEN CHANGED TO HANDLE SCREENED EMISSION YET!!!!
 *
 * rp: Object representing the particle emitting state.
 */
void ConeBremsstrahlungScreenedEmission::CalculateSpectrum(RadiationParticle *rp) {
    slibreal_t
        gamma = rp->GetGamma(),
        gamma2 = rp->GetGamma() * rp->GetGamma(),
        p02 = rp->GetP2(),
        p0 = sqrt(rp->GetP2());
    
    slibreal_t
        pf = this->Z2 * this->r02Alpha,
        TT1 = 4.0 / 3.0,
        TT2, TT3, TT4, TT5,
        L, LT1, LT2, LT3, LT4, LT5, LTF,
        E, E2, p, p2, k, eps, eps0;

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
        eps0 = 2.0*log(gamma+p);

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
 * Integrates the synchrotron spectrum to produce a
 * total emitted power in a given spectral range.
 */
void ConeBremsstrahlungScreenedEmission::IntegrateSpectrum() {
    unsigned int i;
    slibreal_t s = 0.5*(I[0] + I[nwavelengths-1]);

    for (i = 1; i < nwavelengths-1; i++)
        s += I[i];

    this->power = s * (wavelengths[1]-wavelengths[0]);
}

