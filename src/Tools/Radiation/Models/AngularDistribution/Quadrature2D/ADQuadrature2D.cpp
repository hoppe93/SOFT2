/**
 * Implementation of the base class 'ADQuadrature2D'.
 */

#include "Tools/Radiation/Models/AngularDistribution/ADQuadrature2D.h"

using namespace __Radiation;

/**
 * Constructor.
 *
 * a: Emission object that generates the integrand.
 * d: Detector object used in evaluating certain
 *    integration variables.
 */
ADQuadrature2D::ADQuadrature2D(ADEmission *a, Detector *d) {
    this->emission = a;
    this->detector = d;
    this->nwavelengths = d->GetNWavelengths();

    if (nwavelengths > 0) {
        this->I = new slibreal_t[nwavelengths];
        this->Q = new slibreal_t[nwavelengths];
        this->U = new slibreal_t[nwavelengths];
        this->V = new slibreal_t[nwavelengths];
    } else {
        this->I = this->Q = this->U = this->V = nullptr;
    }
}

/**
 * Destructor.
 */
ADQuadrature2D::~ADQuadrature2D() {
    if (this->I != nullptr)
        delete [] this->I;
    if (this->Q != nullptr)
        delete [] this->Q;
    if (this->U != nullptr)
        delete [] this->U;
    if (this->V != nullptr)
        delete [] this->V;
}

/**
 * Evaluate the sines and cosines of the angles
 * 'psi' and 'mu'.
 *
 * dX, dY: Displacement from center of detector of the
 *         point considered.
 * rp:     Particle emitting state.
 * a:      Contains the calculated values of the angles
 *         on return.
 */
void ADQuadrature2D::EvaluateAngles(
    slibreal_t dX, slibreal_t dY,
    RadiationParticle *rp, struct angles &a
) {
    Vector<3> &rcp0 = rp->GetRCP();
    Vector<3> &pHat = rp->GetPHat();
    Vector<3> &e1   = this->detector->GetEHat1();
    Vector<3> &e2   = this->detector->GetEHat2();

    a.n     = rcp0 + dX*e1 + dY*e2;
    a.rcp2  = a.n.Dot(a.n);
    a.n    /= sqrt(a.rcp2);
    a.cosMu =-pHat.Dot(a.n);

    a.cos2Mu = a.cosMu*a.cosMu;
    a.sin2Mu = 1.0 - a.cos2Mu;
    a.sinMu  = sqrt(a.sin2Mu);
}

