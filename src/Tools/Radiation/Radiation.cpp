/**
 * Implementation of the 'Radiation' tool.
 *
 * This tool calculates various forms of radiation from
 * particle and guiding-center orbits.
 */

#include <omp.h>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "SOFT.h"
#include "Tools/Radiation/Radiation.h"

using namespace __Radiation;

/**
 * Constructor.
 */
Radiation::Radiation(
    MagneticField2D *mf, ParticleGenerator *pgen,
    ParticlePusher *pusher
) : Tool("Radiation") {
    this->magfield = mf;
    this->pgen = pgen;
    pusher->ToggleJacobianCalculation(true);
}

/**
 * Destructor.
 */
Radiation::~Radiation() {
    delete [] cosphi;
    delete [] sinphi;
}

/**
 * Apply this tool to the given orbit object.
 *
 * o: Particle/guiding-center orbit to handle.
 * p: Particle object associcated with the orbit.
 */
void Radiation::Handle(Orbit *o, Particle *p) {
    if (this->quadrature == QUADRATURE_FINDSOV)
        HandleTrapzImproved(o, p);
    else
        HandleTrapz(o, p);
}

/**
 * Initialize this tool before using it.
 */
void Radiation::Initialize() {
    for (unsigned int i = 0; i < this->noutput; i++) {
        this->output[i]->Initialize();
    }
}

/**
 * Check if the given point is within the field-of-view.
 *
 * x, y, z: Cartesian coordinates of point to check.
 * rcp:     Vector to store 'r_cp' components in. Contains
 *          vector from detector (camera) to particle on
 *          return.
 */
bool Radiation::IsWithinFieldOfView(slibreal_t x, slibreal_t y, slibreal_t z, Vector<3>& rcp) {
    Vector<3> X = detector->GetPosition();
    Vector<3> n = detector->GetDirection();

    rcp[0] = x - X[0];
    rcp[1] = y - X[1];
    rcp[2] = z - X[2];

    slibreal_t ndotr = n.Dot(rcp);
    slibreal_t fov = rcp.Norm()*detector->GetCosVisionAngleFOV();

    //return (n.Dot(rcp) >= rcp.Norm()*detector->GetCosVisionAngleFOV());
    return (ndotr >= fov);
}

/**
 * Apply all output handlers to the current model state.
 */
void Radiation::RegisterOutput(RadiationParticle *rp) {
    for (unsigned int i = 0; i < noutput; i++) {
        this->output[i]->Handle(this->detector, this->model, rp);
    }
}

