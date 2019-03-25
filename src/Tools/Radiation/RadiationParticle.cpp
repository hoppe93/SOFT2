/**
 * Implementation of the 'RadiationParticle' class.
 */

#include <softlib/config.h>
#include <softlib/Vector.h>
#include "Tools/Radiation/RadiationParticle.h"

using namespace __Radiation;

/**
 * Constructor.
 *
 * o:      Orbit object to initialize the particle from.
 * ti:     Time index at which to evaluate the orbit.
 * detpos: 3-vector denoting position of detector.
 */
RadiationParticle::RadiationParticle(Orbit *o, unsigned int ti, Vector<3>& detpos) :
    RadiationParticle(
        o->GetX(ti), o->GetP(ti), o->GetJdtdrho(ti),
        o->GetPpar(ti), o->GetPperp(ti),
        o->GetGamma(ti), o->GetP2(ti), detpos,
        o->GetBabs(ti), o->GetB(ti), o->GetBhat(ti),
        o->GetMass(), o->GetCharge(), o->GetIndexR(),
        o->GetIndexP1(), o->GetIndexP2(),
        o->GetGradB(ti), o->GetCurlB(ti), o->GetBJacobian(ti)
    ) {}
RadiationParticle::RadiationParticle(
    const Vector<3>& x, const Vector<3>& p, slibreal_t Jdtdrho,
    slibreal_t ppar, slibreal_t pperp,
    slibreal_t gamma, slibreal_t p2, const Vector<3>& detpos,
    slibreal_t B, const Vector<3>& Bvec, const Vector<3>& bHat,
    slibreal_t m, slibreal_t q, unsigned int ir, unsigned int ip1,
    unsigned int ip2, slibreal_t *gradB,
    slibreal_t *curlB, slibreal_t **jacobianB
) {
    this->x = x;
    this->p = p;
    this->Jdtdrho = Jdtdrho;
    this->ppar = ppar;
    this->pperp = pperp;
    this->gamma = gamma;
    this->p2 = p2;
    this->B = B;
    this->m = m;    // Particle mass
    this->q = q;    // Particle charge

    phat = p / p.Norm();
    this->R = hypot(x[0], x[1]);

    this->Bvec = Bvec;
    this->bHat = bHat;
    this->rcp = (x-detpos);
    this->rcp2 = rcp.Dot(this->rcp);
    this->rcplen = sqrt(this->rcp2);
    this->rcphat = this->rcp / this->rcplen;

    this->ir = ir;
    this->ip1 = ip1;
    this->ip2 = ip2;

    if (gradB != nullptr)
        this->gradB = gradB;
    if (curlB != nullptr) {
        this->curlB = curlB;
        if (gradB != nullptr)
            this->curlBhat = (this->curlB + Vector<3>::Cross(this->bHat, this->gradB, this->curlBhat)) / B;
    }
    if (jacobianB != nullptr) {
        this->jacobianB = new slibreal_t*[3];
        this->jacobianB[0] = new slibreal_t[3*3];
        this->jacobianB[1] = this->jacobianB[0]+3;
        this->jacobianB[2] = this->jacobianB[1]+3;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                this->jacobianB[i][j] = jacobianB[i][j];
    }
}

/**
 * Destructor.
 */
RadiationParticle::~RadiationParticle() {
    if (this->jacobianB != nullptr) {
        delete [] this->jacobianB[0];
        delete [] this->jacobianB;
    }
}

/**
 * Sets the differential elements for this particle.
 *
 * dphi: Toroidal step length (angular difference).
 * drho: Step length of radial coordinate.
 * dp1:  Step length of momentum coordinate 1.
 * dp2:  Step length of momentum coordinate 2.
 */
void RadiationParticle::SetDifferentialElements(
    slibreal_t dphi, slibreal_t drho,
    slibreal_t dp1, slibreal_t dp2
) {
    this->dphi = dphi;
    this->Rdphi = this->R*dphi;
    this->drho = drho;
    this->dp1 = dp1;
    this->dp2 = dp2;

    this->diffel = Jdtdrho*dphi * dp1*dp2;
}

/**
 * Update the position and momentum of this particle.
 *
 * x, y:   X and Y coordinates of particle position.
 * px, py: X and Y coordinates of particle momentum.
 * rcp:    Vector from detector to particle.
 * sinr:   Sine of rotation angle.
 * cosr:   Cosine of rotation angle.
 */
void RadiationParticle::UpdateXY(
    slibreal_t x, slibreal_t y, slibreal_t px, slibreal_t py,
    Vector<3>& rcp
) {
    this->x[0] = x;
    this->x[1] = y;
    this->p[0] = px;
    this->p[1] = py;
    this->rcp = rcp;

    this->phat = this->p / this->p.Norm();
}

