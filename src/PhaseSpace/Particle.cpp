/**
 * Implementation of the Particle class.
 */

#include <softlib/config.h>
#include <softlib/Vector.h>

#include "PhaseSpace/Particle.h"

using namespace std;

/**
 * Constructors
 */
Particle::Particle() { }

/**
 * Initializes the particle position.
 * Sets the initial radius.
 *
 * position_type: The type of placement that
 *    was done (Particle::POSITION_???).
 * rho:           Particle radial position.
 * drift_shift:   Orbit drift shift.
 */
void Particle::InitializePosition(const int position_type, const slibreal_t rho, const slibreal_t z0, const slibreal_t drho, const slibreal_t drift_shift) {
	this->position_type = position_type;
	this->rho = rho;
    this->z0  = z0;
	this->drift_shift = drift_shift;
    this->drho = drho;

	this->position_initialized = true;
}

/**
 * Nudge this particle, i.e. move it in the outwards
 * radial direction by the given amount.
 *
 * dr: Amount by which to move the particle in r (in meters).
 */
void Particle::Nudge(const slibreal_t dr, enum nudge_direction nd) {
    if (nd == NUDGE_OUTWARDS)
        this->InitializePosition(position_type, rho+dr, z0, drho, drift_shift);
    else
        this->InitializePosition(position_type, rho-dr, z0, drho, drift_shift);
}

/***********
 * GETTERS *
 ***********/
/**
 * Returns the three-momentum corresponding to
 * the given magnetic field direction and the
 * properties of this particle.
 *
 * bhat: Magnetic-field unit vector.
 *
 * TODO: Use the specified gyrophase.
 */
Vector<3> Particle::Get3Momentum(Vector<3>& bhat) {
    Vector<3> p, bhato;

    if (fabs(bhat[0]) > fabs(bhat[1]) && fabs(bhat[2]) > fabs(bhat[0])) {
        bhato[0] = 1.0;
        bhato[1] = 1.0;
        bhato[2] = -(bhat[0] + bhat[1]) / bhat[2];
    } else if (fabs(bhat[1]) > fabs(bhat[0])) {
        bhato[0] = 1.0;
        bhato[1] = -(bhat[2] + bhat[0]) / bhat[1];
        bhato[2] = 1.0;
    } else {
        bhato[0] = -(bhat[1] + bhat[2]) / bhat[0];
        bhato[1] = 1.0;
        bhato[2] = 1.0;
    }

    bhato.Normalize();
    p = ppar*bhat + pperp*bhato;
    return p;
}

