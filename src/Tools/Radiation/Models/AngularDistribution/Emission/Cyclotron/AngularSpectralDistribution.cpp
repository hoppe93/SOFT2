/**
 * Implementation of spectral distribution
 * of cyclotron radiation
 */

#include "Tools/Radiation/Models/AngularDistribution/ADCyclotronEmission.h"
#include "Tools/Radiation/cyclotron_func.h"
#include <iostream>
#include <algorithm>

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

void ADCyclotronEmission::PreparePolarization(RadiationParticle *rp) {
	PrepareSpectrum(rp);

	this->Efield.nE = this->nwavelengths;

	this->Efield.Ex2  = new slibreal_t[this->Efield.nE];
	this->Efield.Ey2  = new slibreal_t[this->Efield.nE];
	this->Efield.ExEy = new slibreal_t[this->Efield.nE];
}

/**
 * Calculate the minimum and maximum harmonics which emit
 * radiation into the spectral range of the detector.
 */
void ADCyclotronEmission::GetHarmonicsInSpectralRange(
	const slibreal_t cosMu, unsigned int &m_min, unsigned int &m_max
) {
    slibreal_t
		f_max=wavelengths[nwavelengths-1],
		f_min=wavelengths[0];

	// First, determine which harmonics end up in our detector's
	// spectral range.
	
	// rounds down to integer
    m_max=(int)floor(f_max*2*M_PI/omega_0*(1-betapar*cosMu));
	//rounds up to integer
    m_min=(int)ceil(f_min*2*M_PI/omega_0*(1-betapar*cosMu));
	if (m_min == 0)
		m_min = 1;
}

/*
 * Calculate the cyclotron spectrum.
 * The configuration file expects frequency values for @detector -> "spectrum"
 */
void ADCyclotronEmission::CalculateSpectrum(
    Vector<3>&, slibreal_t sinMu,  slibreal_t cosMu
) {
	slibreal_t
		f_max=wavelengths[nwavelengths-1],
		f_min=wavelengths[0],
		prefac  = (cosMu-betapar)/sinMu,
		prefac2 = prefac*prefac,
		// x/w
		x_w=betaperp*sinMu/omega_0,
		// w/m
		w_m=omega_0/(1-betapar*cosMu),
		//spacing between frequencies as stated in the config script
		freq_spacing=(f_max-f_min)/nwavelengths;

    unsigned int m_min=0, m_max=0;
	GetHarmonicsInSpectralRange(cosMu, m_min, m_max);

    for (unsigned int i = 0; i < nwavelengths; i++)
    	I[i] = 0;

	// Iterate over harmonics
    /*for (unsigned int m = m_min; m < m_max+1; ++m) {
    	if (m==0 ||
			((harmonics!=nullptr) &&
			 (find(harmonics, harmonics+harmonics_number, m) == harmonics+harmonics_number)))
			 continue;
	*/
	for (unsigned int m : harmonics) {
		// Harmonic out-of-spectral-range?
		if (m < m_min || m > m_max)
			continue;

    	slibreal_t
			w = m*w_m,
			J = cyclotron_func1(m,w*x_w),
			J_prime = 0.5*(cyclotron_func1(m-1,w*x_w)-cyclotron_func1(m+1,w*x_w));

		// Corresponding index in the frequency range by rounding down to
		// nearest integer
		unsigned int f_index=(int) ((w/(2*M_PI)-f_min)/freq_spacing);

    	I[f_index] = prefactor*w*w*(prefac2*J*J+betaperp*betaperp*J_prime*J_prime);
	}
}

void ADCyclotronEmission::CalculateAngularDistribution(Vector<3>&, slibreal_t, slibreal_t) {
	// TODO
}
void ADCyclotronEmission::CalculatePolarization(
	RadiationParticle *rp, Vector<3>&,
	slibreal_t sinMu,  slibreal_t cosMu
) {
	slibreal_t
		f_max=wavelengths[nwavelengths-1],
		f_min=wavelengths[0],
		prefac = (cosMu-betapar)/sinMu,
		prefac2 = prefac*prefac,
		// x/w
		x_w=betaperp*sinMu/omega_0,
		// w/m
		w_m=omega_0/(1-betapar*cosMu),
		//spacing between frequencies as stated in the config script
		freq_spacing=(f_max-f_min)/nwavelengths;

    unsigned int m_min=0, m_max=0;
	GetHarmonicsInSpectralRange(cosMu, m_min, m_max);

	for (unsigned int i = 0; i < this->nwavelengths; i++) {
		this->I[i] = 0; this->Q[i] = 0;
		this->U[i] = 0; this->V[i] = 0;
		this->Efield.Ex2[i] = 0; this->Efield.Ey2[i] = 0;
		this->Efield.ExEy[i] = 0;
	}

	// Calculate geometry
	Vector<3> bhat = this->magfield->Eval(rp->GetPosition());
	bhat.Normalize();

	slibreal_t r = rp->GetRCP().Norm();
	Vector<3> rcphat = rp->GetRCP() / r;
	slibreal_t nDotB = rcphat.Dot(bhat);

	this->Efield.zhat = rcphat;
	this->Efield.yhat = -Vector<3>::Cross(rcphat, bhat) / (slibreal_t)sqrt(1.0 - nDotB*nDotB);
	this->Efield.xhat = Vector<3>::Cross(this->Efield.yhat, this->Efield.zhat);

	this->Efield.nE = this->nwavelengths;

	// Iterate over harmonics
	unsigned int prev_f=0;
	for (unsigned int m : harmonics) {
		// Harmonic out-of-spectral-range?
		if (m < m_min || m > m_max)
			continue;

    	slibreal_t
			w = m*w_m,
			J = cyclotron_func1(m,w*x_w),
			J_prime = 0.5*(cyclotron_func1(m-1,w*x_w)-cyclotron_func1(m+1,w*x_w));

		Efield.Ex2[prev_f] = 0;
		Efield.Ey2[prev_f] = 0;
		Efield.ExEy[prev_f] = 0;

		// Corresponding index in the frequency range by rounding down to
		// nearest integer
		unsigned int f_index=(int) ((w/(2*M_PI)-f_min)/freq_spacing);

		Efield.Ex2[f_index]  = prefactor*w*w * prefac2*J*J;
		Efield.Ey2[f_index]  = prefactor*w*w * betaperp*betaperp*J_prime*J_prime;
		Efield.ExEy[f_index] = prefactor*w*w * prefac*betaperp*J*J_prime;

		this->detector->GetOptics()->ApplyOptics(Efield, I, Q, U, V);
	}
}
void ADCyclotronEmission::InitializeToroidalStep(const slibreal_t, const slibreal_t){
	// Nothing to do...
}

