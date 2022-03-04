/**
 * General interface for calculation of CYCLOTRON emission
 */

#include "SOFT.h"
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/AngularDistribution/ADCyclotronEmission.h"

using namespace std;
using namespace __Radiation;

/**
 * Constructor.
 *
 * det:        Detector object defining the detector properties.
 * mf:         Magnetic field object.
 * globset:    Global settings object (used to determine whether
 *             or not to enable drifts; drifts are disabled at the moment).
 */
ADCyclotronEmission::ADCyclotronEmission(
	Detector *det, MagneticField2D *mf, struct global_settings*, 
	int *harmonics_list, int harmonics_no
) : ADEmission(det, mf) {

	if (det->GetNWavelengths() == 0)
		throw ADCyclotronException("Wavelength range needs to be specified?");

	for (int i = 0; i < harmonics_no; i++)
		this->harmonics.push_back(harmonics_list[i]);

	std::sort(this->harmonics.begin(), this->harmonics.end());
}

/**
 * Destructor.
 */
ADCyclotronEmission::~ADCyclotronEmission() { }

/**
 * High-level function for calculating the spectral distribution
 * of Cyclotron radiation from the given particle (state).
 *
 * sinPsi, cosPsi: Angle between _particle_ velocity and observation direction.
 * sinMu,  cosMu:  Angle between _GC_ velocity and observation direction.
 * pol: polarization: not regarded here
 */
slibreal_t ADCyclotronEmission::Evaluate(
	RadiationParticle *rp, Vector<3> &n,slibreal_t sinMu,
	slibreal_t cosMu, bool pol
) {
	if (pol) {
		CalculatePolarization(rp, n, sinMu, cosMu);
		IntegrateSpectrumStokes();
		IntegrateSpectrum();
	} else {
		CalculateSpectrum(n, sinMu, cosMu);
		IntegrateSpectrum();
	}

    return this->power;
}

/**
 * Prepare to handle a particle. Calculates necessary
 * cache quantities.
 *
 * rp:  Particle emission state.
 * pol: polarization: not regarded here
 */
void ADCyclotronEmission::Prepare(RadiationParticle *rp, bool pol) {
	if (pol)
		PreparePolarization(rp);
	else
		PrepareSpectrum(rp);
}

/**
 * Integrate the calculated spectrum.
 */
void ADCyclotronEmission::IntegrateSpectrum() {
    unsigned int i;
    slibreal_t s = 0.5*(I[0] + I[nwavelengths-1]);

    for (i = 1; i < nwavelengths-1; i++)
        s += I[i];

    this->power = s * (wavelengths[1]-wavelengths[0]);
}

/**
 * Integrate Stokes parameters.
 */
void ADCyclotronEmission::IntegrateSpectrumStokes() {
    unsigned x;
    slibreal_t
        q=0.5 * (Q[0] + Q[nwavelengths-1]),
        u=0.5 * (U[0] + U[nwavelengths-1]),
        v=0.5 * (V[0] + V[nwavelengths-1]);

    for (x = 1; x < nwavelengths-1; x++)
        q += Q[x];
    for (x = 1; x < nwavelengths-1; x++)
        u += U[x];
    for (x = 1; x < nwavelengths-1; x++)
        v += V[x];

    this->powerQ = q;
    this->powerU = u;
    this->powerV = v;
}

