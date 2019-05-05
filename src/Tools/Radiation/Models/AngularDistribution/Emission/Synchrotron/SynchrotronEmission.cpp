/**
 * General interface for calculation of synchrotron
 * emission in the angular distribution model.
 */

#include "SOFT.h"
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"

using namespace std;
using namespace __Radiation;

/**
 * Constructor.
 *
 * det:        Detector object defining the detector properties.
 * mf:         Magnetic field object.
 * globset:    Global settings object (used to determine whether
 *             or not to enable drifts).
 * I, Q, U, V: Pre-allocated arrays consisting of 'nwavelengths'
 *             number of elements each, where 'nwavelengths' is
 *             the number of wavelengths to sample in the
 *             spectrum (as set in 'det').
 */
ADSynchrotronEmission::ADSynchrotronEmission(
    Detector *det, MagneticField2D *mf, struct global_settings *globset,
    size_t qagsLimit, slibreal_t qagsEpsRel
) : ADEmission(det, mf) {
    this->includeDrifts = globset->include_drifts;

    if (includeDrifts && det->GetNWavelengths() != 0)
        throw ADSynchrotronException("The model for angular and spectral distribution of synchrotron radiation has no support for drifts.");

    if (includeDrifts) {
        this->qagsLimit = qagsLimit;
        this->qagsEpsRel = (double)qagsEpsRel;
        this->qagsWorkspace = gsl_integration_workspace_alloc(qagsLimit);
    }
}

/**
 * Destructor.
 */
ADSynchrotronEmission::~ADSynchrotronEmission() { }

/**
 * High-level function for calculating the angular distribution
 * of synchrotron radiation from the given particle (state).
 * 
 * sinPsi, cosPsi: Angle between _particle_ velocity and observation direction.
 * sinMu,  cosMu:  Angle between _GC_ velocity and observation direction.
 * pol: If true, calculates the polarization components of
 *      synchrotron radiation.
 */
slibreal_t ADSynchrotronEmission::Evaluate(
    RadiationParticle *rp, Vector<3> &n,
    slibreal_t sinMu,  slibreal_t cosMu, bool pol
) {
    if (nwavelengths == 0) {
        CalculateAngularDistribution(n, sinMu, cosMu);
    } else {
        if (!pol)
            CalculateSpectrum(n, sinMu, cosMu);
        else {
            CalculatePolarization(rp, n, sinMu, cosMu);
            IntegrateSpectrumStokes();
        }

        IntegrateSpectrum();
    }

    return this->power;
}

/**
 * Initialize the current toroidal step.
 */
void ADSynchrotronEmission::InitializeToroidalStep(const slibreal_t sinphi, const slibreal_t cosphi) {
    if (nwavelengths == 0)
        InitializeToroidalStepAD(sinphi, cosphi);
}

/**
 * Prepare to handle a particle. Calculates necessary
 * cache quantities.
 *
 * rp:  Particle emission state.
 * pol: If true, calculates the polarization components of
 *      synchrotron radiation.
 */
void ADSynchrotronEmission::Prepare(RadiationParticle *rp, bool pol) {
    if (nwavelengths == 0)
        PrepareAngularDistribution(rp);
    else {
        if (!pol)
            PrepareSpectrum(rp);
        else
            PreparePolarization(rp);
    }
}

/**
 * Integrate the calculated spectrum.
 */
void ADSynchrotronEmission::IntegrateSpectrum() {
    unsigned int i;
    slibreal_t s = 0.5*(I[0] + I[nwavelengths-1]);

    for (i = 1; i < nwavelengths-1; i++)
        s += I[i];

    this->power = s * (wavelengths[1]-wavelengths[0]);
}

/**
 * Integrate Stokes parameters.
 */
void ADSynchrotronEmission::IntegrateSpectrumStokes() {
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

