/**
 * Module for evaluating screened bremsstrahlung emission in
 * the cone model.
 */

#include "Tools/Radiation/Models/Cone/ConeBremsstrahlungScreenedEmission.h"
#include "Tools/Radiation/RadiationParticle.h"

using namespace __Radiation;
ConeBremsstrahlungScreenedEmission::ConeBremsstrahlungScreenedEmission(Detector *det, MagneticField2D *mf, unsigned int nspecies, slibreal_t *Z, slibreal_t *Z0, slibreal_t *density)
                : ConeEmission(det, mf), nspecies(nspecies)
{
    this->Z = new slibreal_t[nspecies]; 
    this->Z0 = new slibreal_t[nspecies];
    this->density = new slibreal_t[nspecies];
    for(unsigned int i = 0; i < nspecies; i++){
        this->Z[i] = Z[i];
        this->Z0[i] = Z0[i];
        this->density[i] =  density[i];
    }
    this->r02 = ELECTRON_CHARGE*ELECTRON_CHARGE*ELECTRON_CHARGE*ELECTRON_CHARGE/(16*M_PI*M_PI*EPS0*EPS0*LIGHTSPEED*LIGHTSPEED*LIGHTSPEED*LIGHTSPEED*ELECTRON_MASS*ELECTRON_MASS);
}



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
    //printf("%u", nspecies);
    for(unsigned int i=0; i < nspecies; i++){
        I = I + density[i]*Calculate4BS(Z[i]);
        //printf("n = %e, Z = %e \n", density[i], Z[i]);
    }
    this->power = I;
}

/**
 * Calculate eq. 4BS (Koch & Motz) for a given species
 */

slibreal_t ConeBremsstrahlungScreenedEmission::Calculate4BS(slibreal_t Z) {
    slibreal_t Z2fakt = 4*Z*Z*r02/137;
    slibreal_t lnfakt = log(183) - 0.5*log(Z) + 0.0555555555555556;
    
    return Z2fakt*lnfakt;
}

/**
 * Calculate the bremsstrahlung screned emissiion spectrum.
 *
 * rp: Object representing the particle emitting state.
 */
/*
void ConeBremsstrahlungScreenedEmission::CalculateSpectrum(RadiationParticle *rp) {
    slibreal_t
        gamma = rp->GetGamma(),
        gamma2 = rp->GetGamma() * rp->GetGamma()
        //p2 = rp->GetP2(),
        p = sqrt(rp->GetP2());

    slibreal_t m = ELECTRON_MASS,
        c = LIGHTSPEED;        
        
    slibreal_t spec_cont = 0,
        Spec_Int1,
        Spec_Int2;
	
    
    slibreal_t I[nwavelengths] = {0};
    slibreal_t q0;
    slibreal_t k_normed;
    slibreal_t gmkn;
    slibreal_t = pre_fact = 4*r02*alpha;
    
    for (unsigned int i = 0; i < nwavelengths; i++){
        k_normed = 2*M_PI*HBAR/(wavelength[i]*m*c); //factor 1/(m*c) for normalization
        gmkn = gamma-k_normed;
        q0 = p - sqrt(gmkn*gmkn-1)-k_normed;
        
        for (unsigned int j = 0; j < nspecies; j++){
            Spec_Int1 = FirstSpectrumIntegral(Z[j], Z0[j], q0); 
            Spec_Int2 = SecondSpectrumIntegral(Z[j], Z0[j], q0); 
            spec_cont = density[j]*((1+gmkn*gmkn/gamma2)*Spec_Int1 - 2*gmkn/(3*gamma)*Spec_Int2);
        }
        I[i] = pre_fact * wavelength[i]/(2*M_PI*hbar)*spec_cont;
        spec_cont = 0;
        
    }
   
}

slibreal_t ConeBremsstrahlungScreenedEmission::FirstSpectrumIntegral(slibreal_t Z, slibreal_t Z0, slibreal_t q0){
    //kalla på CalculateFormFactor, integrera
    slibreal_t Integral = 0;
    slibreal_t Fjq = CalculateFormFactor(Z, Z0, q);
    return Z*Z + Integral;
}

slibreal_t ConeBremsstrahlungScreenedEmission::SecondSpectrumIntegral(slibreal_t Z, slibreal_t Z0, slibreal_t q0){
    //kalla på CalculateFormFactor, integrera
    slibreal_t Integral = 0;
    slibreal_t Fjq = CalculateFormFactor(Z, Z0, q);
    return Z*Z + Integral;
}
*/

/*
 * Calculates the formfactor Fj(q), only very basic fomrula so far
*/

slibreal_t ConeBremsstrahlungScreenedEmission::CalculateFormFactor(slibreal_t Z, slibreal_t Z0, slibreal_t q){
    slibreal_t Nej2 = (Z - Z0)*(Z - Z0);
    slibreal_t a_bar = pow(9*M_PI*Nej2, 0.3333333333333333333333333)/(2*alpha*Z); //Kirilov, from "generalized collision..."
    slibreal_t qabar3o2 = sqrt(q*q*q*a_bar*a_bar*a_bar);
    return Nej2/(1+qabar3o2);
    
}


/**
 * Calculate the bremsstrahlung spectrum. OLD!!!! NOT SCREENED EMISSION!
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

