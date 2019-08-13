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
    qagsWS = gsl_integration_workspace_alloc(qagsLimit); //GSL workspace
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
 */
void ConeBremsstrahlungScreenedEmission::CalculateTotalEmission() {
    slibreal_t I = 0;   
    for(unsigned int i=0; i < nspecies; i++){
        I = I + density[i]*Calculate4BS(Z[i]);
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

void ConeBremsstrahlungScreenedEmission::CalculateSpectrum(RadiationParticle *rp) {
    slibreal_t
        gamma = rp->GetGamma(),
        gamma2 = rp->GetGamma() * rp->GetGamma(),
        p = sqrt(rp->GetP2());
        
    slibreal_t spec_cont = 0,
        Spec_Int1 = 0,
        Spec_Int2 = 0,
        Nej,
        a_bar;
    slibreal_t q0,
        k_normed,
        gmkn,
        pre_fact = 4*r02*alpha;
    //printf("c");
    for(unsigned int i = 0; i < nwavelengths; i++){
       I[i] = 0;
    }
    for (unsigned int j = 0; j < nspecies; j++){
            Nej = (Z[j] - Z0[j]);
            a_bar = pow(9*M_PI*Nej*Nej, 0.3333333333333333333333333)/(2*alpha*Z[j]);
    
		for (unsigned int i = 0; i < nwavelengths; i++){
		    k_normed = wavelengths[i]; 
		    if (k_normed > gamma - 1) { //VAD ska det det vara?
		        continue;
		    }

		    gmkn = gamma-k_normed;
		    q0 = p - sqrt(gmkn*gmkn-1)-k_normed;
			
		   
		    Spec_Int1 = FirstSpectrumIntegral(Z[j], Nej, a_bar, q0); 
		    Spec_Int2 = SecondSpectrumIntegral(Z[j], Nej, a_bar, q0); 
		    spec_cont = density[j]*((1+gmkn*gmkn/gamma2)*(Z[j]*Z[j] + Spec_Int1) - 2*gmkn/(3*gamma)*(5*Z[j]*Z[j]/6 + Spec_Int2));
		    
		    I[i] = I[i] + pre_fact*spec_cont/wavelengths[i];
		     
		} 
    }
}

//First spectrum integral
slibreal_t ConeBremsstrahlungScreenedEmission::FirstSpectrumIntegral(slibreal_t Z, slibreal_t Nej, slibreal_t a_bar, slibreal_t q0){
    
    struct func_params params = {Z, Nej, a_bar, q0, 0};
    gsl_function F;
    F.function = &First_Integrand;
    F.params = (void*)&params;
    slibreal_t result,
        qagsAbsErr;

    gsl_integration_qags(&F, q0, 1, this->qagsEpsAbs, this->qagsEpsRel, this->qagsLimit, qagsWS, &result, &qagsAbsErr);
    
    return result;
}

//second spectrum integral
slibreal_t ConeBremsstrahlungScreenedEmission::SecondSpectrumIntegral(slibreal_t Z, slibreal_t Nej, slibreal_t a_bar, slibreal_t q0){
   
    struct func_params p = {Z, Nej, a_bar, q0, 0}; //Function Parameters
    gsl_function F;
    
    if(q0 <= 0.05){
		F.function = &Second_Integrand_app_sq0;
	}
	else{
		p = {Z, Nej, a_bar, q0, log(q0)};
    	F.function = &Second_Integrand_app_lq0;
	}
	
	//F.function = &Second_Integrand; //For no approximation in second integrand
    F.params = (void*)&p;
    slibreal_t result,
        qagsAbsErr;
    
    gsl_integration_qags(&F, q0, 1, this->qagsEpsAbs, this->qagsEpsRel, this->qagsLimit, qagsWS, &result, &qagsAbsErr);
    return result;
}

//First spectrum integrand
double ConeBremsstrahlungScreenedEmission::First_Integrand(slibreal_t q, void *p){
    struct func_params *params = (struct func_params*)p; 
    slibreal_t Fq = CalculateFormFactor(params->Nej, params->a_bar, q);
    slibreal_t Z = params->Z,
        q0 = params->q0,
        ZmF2 = (Z-Fq)*(Z-Fq),
        q_fun = (q-q0)*(q-q0)/(q*q*q);
    
    return ZmF2*q_fun;
}

//second spectrum integrand, no approximation
double ConeBremsstrahlungScreenedEmission::Second_Integrand(slibreal_t q, void *p){
    struct func_params *params = (struct func_params*)p;
    slibreal_t Fq = CalculateFormFactor(params->Nej, params->a_bar, q),
        q0 = params->q0,
        Z = params->Z,
        ZmF2 = (Z-Fq)*(Z-Fq),
        q3 = q*q*q,
        q02 = q0*q0,
        q_fun = (q3 + 3*q*q02*(1 - 2*q*q02*log(q/q0))-4*q0*q02)/(q3*q);
    
    return ZmF2*q_fun;
}

//Second spectrum integrand, approximation of log-term for q0 <= 0.05
double ConeBremsstrahlungScreenedEmission::Second_Integrand_app_sq0(slibreal_t q, void *p){
    struct func_params *params = (struct func_params*)p;
    slibreal_t Fq = CalculateFormFactor(params->Nej, params->a_bar, q),
        q0 = params->q0,
        Z = params->Z,
        ZmF2 = (Z-Fq)*(Z-Fq),
        q3 = q*q*q,
        q02 = q0*q0,
        q_fun = (q3 + 3*q*q02 - 4*q0*q02)/(q3*q);
    
    return ZmF2*q_fun;
}

//Second spectrum integrand, approximation of log-term for q0 > 0.05
double ConeBremsstrahlungScreenedEmission::Second_Integrand_app_lq0(slibreal_t q, void *p){
    struct func_params *params = (struct func_params*)p;
		//printf("lnq0 = %e \n", params->lnq0);
    slibreal_t Fq = CalculateFormFactor(params->Nej, params->a_bar, q),
        q0 = params->q0,
        Z = params->Z,
        ZmF2 = (Z-Fq)*(Z-Fq),
        q3 = q*q*q,
        q02 = q0*q0,
		qm1 = q-1, qm2 = qm1*qm1,
        q_fun = (q3 + 3*q*q02*(1 - 2*q*q02*(qm1 - qm2/2 - params->lnq0))-4*q0*q02)/(q3*q);
    	//ln(q) is here approximated to second order
    return ZmF2*q_fun;
}

//Calculates the formfactor Fj(q), analythical formula from Kirilov

slibreal_t ConeBremsstrahlungScreenedEmission::CalculateFormFactor(slibreal_t Nej, slibreal_t a_bar, slibreal_t q){
    slibreal_t qabar3o2 = sqrt(q*q*q*a_bar*a_bar*a_bar);
    return Nej/(1+qabar3o2); 
}

/**
 * Integrates the synchrotron spectrum to produce a
 * total emitted power in a given spectral range.
 */
void ConeBremsstrahlungScreenedEmission::IntegrateSpectrum() {
    unsigned int i = 1;
    slibreal_t s = 0.5*(I[0] + I[nwavelengths-1]);
	//printf("nwl = %u \n I = %e \n", nwavelengths, I[0]);
    for (i = 1; i < nwavelengths-1; i++){
        s += I[i];
		//printf("int_i = %u \n wl = %e \n I = %e \n", i, wavelengths[i], I[i]);
	}

    this->power = s * (wavelengths[1]-wavelengths[0]);
}

