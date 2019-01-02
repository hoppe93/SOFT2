/**
 * Implementation of the Orbit class.
 */

#include <cstring>
#include <softlib/config.h>
#include "Orbit/Orbit.h"

#include <omp.h>
#include "Orbit/GuidingCenterEquation.h"

/**
 * Constructor.
 *
 * nt: Number of time points to allocate.
 */
Orbit::Orbit(unsigned int nt, bool calcBDerivatives) {
    this->ntau = nt;
    this->hasBDerivatives = calcBDerivatives;

    this->tau        = new slibreal_t[nt];
    this->x          = new slibreal_t[nt*3];
    this->p          = new slibreal_t[nt*3];
    this->ppar       = new slibreal_t[nt];
    this->pperp      = new slibreal_t[nt];
    this->Jdtdrho    = new slibreal_t[nt];
    this->Jp         = new slibreal_t[nt];
    this->_solution  = new slibreal_t[nt*6];
    this->_solution2 = new slibreal_t[nt*6];

    this->Babs    = new slibreal_t[nt];
    this->B       = new slibreal_t[nt*3];
    this->Beffpar = new slibreal_t[nt];
    this->bhat    = new slibreal_t[nt*3];
    this->p2      = new slibreal_t[nt];
    this->ppar2   = new slibreal_t[nt];
    this->pperp2  = new slibreal_t[nt];
    this->gamma   = new slibreal_t[nt];
    
    if (calcBDerivatives) {
        this->gradB = new slibreal_t[nt*3];
        this->curlB = new slibreal_t[nt*3];
        this->jacobianB = new slibreal_t*[nt*3];
        this->jacobianB[0] = new slibreal_t[nt*3*3];

        for (unsigned int i = 1; i < nt*3; i++)
            this->jacobianB[i] = this->jacobianB[i-1]+3;
    } else {
        this->gradB = nullptr;
        this->curlB = nullptr;
        this->jacobianB = nullptr;
    }
}

/**
 * Destructor.
 */
Orbit::~Orbit() {
    if (this->hasBDerivatives) {
        delete [] this->jacobianB[0];
        delete [] this->jacobianB;
        delete [] this->curlB;
        delete [] this->gradB;
    }

    delete [] this->gamma;
    delete [] this->pperp2;
    delete [] this->ppar2;
    delete [] this->p2;
    delete [] this->bhat;
    delete [] this->B;
    delete [] this->Babs;
    delete [] this->Beffpar;

    delete [] this->_solution2;
    delete [] this->_solution;
    delete [] this->Jp;
    delete [] this->Jdtdrho;
    delete [] this->pperp;
    delete [] this->ppar;
    delete [] this->p;
    delete [] this->x;
    delete [] this->tau;

    this->ntau = 0;
}

/**
 * Copies the fields of this Orbit object into
 * the fields of the given orbit object. Note that
 * pointers are _not_ copied, but rather the values
 * they point to.
 *
 * o: Orbit object to copy this object into.
 */
void Orbit::CopyTo(Orbit *o) {
    if (o->GetNTau() != ntau)
        throw OrbitException("Orbit::CopyTo(): Orbit object is not of the correct size.");

    memcpy(o->GetTau(), tau, ntau*sizeof(slibreal_t));
    memcpy(o->GetX(), x, 3*ntau*sizeof(slibreal_t));
    memcpy(o->GetP(), p, 3*ntau*sizeof(slibreal_t));
    memcpy(o->GetPpar(), ppar, ntau*sizeof(slibreal_t));
    memcpy(o->GetPperp(), pperp, ntau*sizeof(slibreal_t));
    memcpy(o->GetJdtdrho(), Jdtdrho, ntau*sizeof(slibreal_t));
    memcpy(o->GetJp(), Jp, ntau*sizeof(slibreal_t));
    memcpy(o->GetInternalSolution(), _solution, 6*ntau*sizeof(slibreal_t));
    memcpy(o->GetInternalSolutionSecondary(), _solution2, 6*ntau*sizeof(slibreal_t));
    memcpy(o->GetBabs(), Babs, ntau*sizeof(slibreal_t));
    memcpy(o->GetB(), B, 3*ntau*sizeof(slibreal_t));
    memcpy(o->GetBeffpar(), Beffpar, ntau*sizeof(slibreal_t));
    memcpy(o->GetBhat(), bhat, 3*ntau*sizeof(slibreal_t));
    memcpy(o->GetP2(), p2, ntau*sizeof(slibreal_t));
    memcpy(o->GetPpar2(), ppar2, ntau*sizeof(slibreal_t));
    memcpy(o->GetPperp2(), pperp2, ntau*sizeof(slibreal_t));
    memcpy(o->GetGamma(), gamma, ntau*sizeof(slibreal_t));

    this->hasBDerivatives = o->HasBDerivatives();
    if (this->hasBDerivatives) {
        memcpy(o->GetGradB(), gradB, 3*ntau*sizeof(slibreal_t));
        memcpy(o->GetCurlB(), curlB, 3*ntau*sizeof(slibreal_t));
        memcpy(o->GetBJacobian(), jacobianB, 3*3*ntau*sizeof(slibreal_t));
    }

    this->ir = o->GetIndexR();
    this->ip1 = o->GetIndexP1();
    this->ip2 = o->GetIndexP2();
}

/**
 * Reads data from the given integrator object and stores
 * the corresponding orbit (obtained through calls to the
 * 'SOFTEquation' object).
 *
 * tend:  Last time point at which to evaluate the orbit.
 * intg1: Integrator object containing solution (6-dimensional).
 * intg2: Second integrator object containing the solution of a
 *        (very) slightly different orbit, that can be used to
 *        calculate the spatial Jacobian determinant.
 * eqn:   Pointer to SOFTEquation solved by the integrator.
 * nudge: Nudge value used to generate 'intg2'.
 */
Orbit *Orbit::Create(slibreal_t tend, Integrator<6> *intg1, Integrator<6> *intg2, SOFTEquation *eqn, Particle *p, slibreal_t nudge, orbit_class_t cl) {
    this->ir = p->GetIndexR();
    this->ip1 = p->GetIndexP1();
    this->ip2 = p->GetIndexP2();
    this->m = p->GetMass();
    this->q = p->GetCharge();

    intg1->OutputDense(
        this->ntau,
        0.0, tend,
        this->_solution, this->tau
    );

    // Evaluate secondary orbit if available
    if (intg2 != nullptr) {
        intg2->OutputDense(
            this->ntau,
            0.0, tend,
            this->_solution2, this->tau
        );

        eqn->ToOrbitQuantities(this->_solution, this->_solution2, this, nudge, cl);
    } else
        eqn->ToOrbitQuantities(this->_solution, nullptr, this, nudge, cl);

    return this;
}

