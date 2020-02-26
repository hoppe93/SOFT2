/**
 * Implementation of spectral distribution
 * of cyclotron radiation
 */

#include "Tools/Radiation/Models/AngularDistribution/ADCyclotronEmission.h"
#include "Tools/Radiation/cyclotron_func.h"
#include <iostream>

using namespace std;
using namespace __Radiation;

/*************************
 * PREPARATION FUNCTIONS *
 *************************/

void ADCyclotronEmission::PrepareSpectrum(RadiationParticle *rp) {
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


    this->beta = sqrt(beta2);
    this->cosThetap = rp->GetPpar() / p;
    this->sinThetap = rp->GetPperp() / p;
    this->betapar = beta*cosThetap;

    this->betaperp= beta*sinThetap;
    this->omega_0=q*B /m /gamma;
    this->prefactor = q*q/(8*M_PI*M_PI*EPS0*LIGHTSPEED);

    this->lambdac = 4.0*M_PI*LIGHTSPEED*gammapar*igamma2*m / (3.0*q*B);
}

/*
 * Calculate the cyclotron spectrum.
 * The configuration file expects frequency values for @detector -> "spectrum"
 */

void ADCyclotronEmission::CalculateSpectrum(
    Vector<3> &n, slibreal_t sinMu,  slibreal_t cosMu
) {

    unsigned int m_min,m_max;
    unsigned int m;
    slibreal_t
				prefac2= pow((cosMu-betapar)/sinMu,2),
				x_w=betaperp*sinMu/omega_0,			// x/w
    			w_m=omega_0/(1-betapar*cosMu),		// w/m
				f_max=wavelengths[nwavelengths-1],
				f_min=wavelengths[0],
    			freq_spacing=(f_max-f_min)/nwavelengths;		//spacing between frequencies as stated in the config script

    m_max=(int) floor(f_max*2*M_PI/omega_0*(1-betapar*cosMu)); // rounds up to integer
    m_min=(int) ceil(f_min*2*M_PI/omega_0*(1-betapar*cosMu)); //rounds up to integer

    for (unsigned int i = 0; i < nwavelengths; i++)
    	I[i] = 0;


    for (m = m_min; m < m_max+1; ++m) {
    	if (m==0)continue;

		slibreal_t
		w= m*w_m,
		J = cyclotron_func1(m,w*x_w),
		J_prime = 0.5*(cyclotron_func1(m-1,w*x_w)-cyclotron_func1(m+1,w*x_w));
		unsigned int f_index=(int) ((w/(2*M_PI)-f_min)/freq_spacing);			//corresponding index in the frequency range by rounding down to nearest int

    	I[f_index]=prefactor*w*w*(prefac2*J*J+betaperp*betaperp*J_prime*J_prime);
	}

}


void ADCyclotronEmission::CalculateAngularDistribution(Vector<3>&, slibreal_t, slibreal_t) {

}
void ADCyclotronEmission::CalculatePolarization(RadiationParticle*, Vector<3>&, slibreal_t, slibreal_t){

}
void ADCyclotronEmission::InitializeToroidalStep(const slibreal_t, const slibreal_t){

}

