/**
 * Implementation of the 'Orbits' tool.
 */

#include <cstring>
#include "Orbit/Orbit.h"
#include "SOFTException.h"
#include "Orbit/ParticlePusher.h"
#include "Tools/Orbits.h"

slibreal_t
    **Orbits::tau,		/* Time points (orbit time) */
    **Orbits::x,		/* Particle/guiding-center position vector */
    **Orbits::p,		/* Particle/guiding-center momentum vector */
    **Orbits::ppar,		/* Parallel momentum of particle */
    **Orbits::pperp,	/* Perpendicular momentum of particle */
    **Orbits::Jdtdrho,	/* Trajectory coordinate Jacobian (times dtau * drho) */
    **Orbits::solution; /* Temporary storage of integrator solution */

/* Related vectors */
slibreal_t
    **Orbits::Babs,      /* Magnetic field strength vector (1-dimensional) */
    **Orbits::B,         /* Magnetic field vector (3-dimensional) */
    **Orbits::bhat,      /* Magnetic field unit vector (3-dimensional) */
    **Orbits::p2,        /* Momentum squared */
    **Orbits::ppar2,     /* Parallel momentum squared */
    **Orbits::pperp2,    /* Perpendicular momentum squared */
    **Orbits::gamma;     /* Lorentz factor (or energy) */

slibreal_t
    *Orbits::drift_shift,/* Orbit drift shift in outer midplane */
	*Orbits::f;			 /* Distribution function */

orbit_class_t *Orbits::classification;  /* Orbit classifications */

/**
 * Constructor.
 *
 * Initialize the Orbits tool and adjust it
 * to suit the configured phase-space.
 */
Orbits::Orbits(MagneticField2D *mf, ParticleGenerator *pg, ParticlePusher *pusher)
        : Tool("Orbits"), OutputModule(mf, pg) {

    nr  = pg->GetNr();
    np1 = pg->GetN1();
    np2 = pg->GetN2();
    norbits = nr*np1*np2;
    ntau = pusher->GetNTimeSteps();
    this->pusher = pusher;
    
    #pragma omp critical (Orbits_Init)
    {
        if (tau == nullptr) {
            this->allocatedBytes = 0;

            Orbits::tau      = Allocate(norbits, ntau, 1);
            Orbits::x        = Allocate(norbits, ntau, 3);
            Orbits::p        = Allocate(norbits, ntau, 3);
            Orbits::ppar     = Allocate(norbits, ntau, 1);
            Orbits::pperp    = Allocate(norbits, ntau, 1);
            Orbits::Jdtdrho  = Allocate(norbits, ntau, 1);
            Orbits::solution = Allocate(norbits, ntau, 6);

            Orbits::Babs     = Allocate(norbits, ntau, 1);
            Orbits::B        = Allocate(norbits, ntau, 3);
            Orbits::bhat     = Allocate(norbits, ntau, 3);
            Orbits::p2       = Allocate(norbits, ntau, 1);
            Orbits::ppar2    = Allocate(norbits, ntau, 1);
            Orbits::pperp2   = Allocate(norbits, ntau, 1);
            Orbits::gamma    = Allocate(norbits, ntau, 1);

            Orbits::classification = Allocate_class(norbits);
            Orbits::drift_shift    = Allocate_prop(norbits);
			Orbits::f              = Allocate_prop(norbits);
        }
    }
}

/**
 * Allocate an array for storing 'norbits'
 * 'dims'-dimensional data, consisting of
 * 'ntau' points each.
 *
 * norbits: Number of orbits to allocate space for (first dimension).
 * ntau:    Number of timesteps to allocate space for.
 * dims:    Number of dimensions in the quantity (i.e., for position,
 *          represented as a 3-vector, this is 3).
 */
slibreal_t **Orbits::Allocate(const unsigned int norbits, const unsigned int ntau, const unsigned int dims) {
    slibreal_t **p = new slibreal_t*[norbits];
    p[0] = new slibreal_t[norbits*ntau*dims];

    this->allocatedBytes += norbits*ntau*dims * sizeof(slibreal_t);

    for (unsigned int i = 1; i < norbits; i++)
        p[i] = p[i-1] + ntau*dims;

    return p;
}

/**
 * Allocate an array for storing 'norbits'
 * orbit classification values.
 *
 * norbits: Number of orbits to allocate space for.
 */
orbit_class_t *Orbits::Allocate_class(const unsigned int norbits) {
    orbit_class_t *p = new orbit_class_t[norbits];
    this->allocatedBytes += norbits * sizeof(orbit_class_t);

    for (unsigned int i = 0; i < norbits; i++) {
        p[i] = ORBIT_CLASS_DISCARDED;
    }

    return p;
}

/**
 * Allocate an array for storing 'norbits'
 * orbit property values (1D arrays).
 *
 * norbits: Number of orbits to allocate space for.
 */
slibreal_t *Orbits::Allocate_prop(const unsigned int norbits) {
    slibreal_t *p = new slibreal_t[norbits];
    this->allocatedBytes += norbits * sizeof(slibreal_t);

    return p;
}

/**
 * Handle the given orbit.
 * 
 * o: Orbit to handle.
 */
void Orbits::Handle(Orbit *o, Particle *part) {
    unsigned int
        ir = o->GetIndexR(),
        ip1 = o->GetIndexP1(),
        ip2 = o->GetIndexP2(),
        index = ir*np1*np2 + ip1*np2 + ip2;

    memcpy(Orbits::tau[index], o->GetTau(), ntau*sizeof(slibreal_t));
    memcpy(Orbits::x[index], o->GetX(), 3*ntau*sizeof(slibreal_t));
    memcpy(Orbits::p[index], o->GetP(), 3*ntau*sizeof(slibreal_t));
    memcpy(Orbits::ppar[index], o->GetPpar(), ntau*sizeof(slibreal_t));
    memcpy(Orbits::pperp[index], o->GetPperp(), ntau*sizeof(slibreal_t));
    memcpy(Orbits::Jdtdrho[index], o->GetJdtdrho(), ntau*sizeof(slibreal_t));
    memcpy(Orbits::solution[index], o->GetInternalSolution(), 6*ntau*sizeof(slibreal_t));

    memcpy(Orbits::Babs[index], o->GetBabs(), ntau*sizeof(slibreal_t));
    memcpy(Orbits::B[index], o->GetB(), 3*ntau*sizeof(slibreal_t));
    memcpy(Orbits::bhat[index], o->GetBhat(), 3*ntau*sizeof(slibreal_t));
    memcpy(Orbits::p2[index], o->GetP2(), ntau*sizeof(slibreal_t));
    memcpy(Orbits::ppar2[index], o->GetPpar2(), ntau*sizeof(slibreal_t));
    memcpy(Orbits::pperp2[index], o->GetPperp2(), ntau*sizeof(slibreal_t));
    memcpy(Orbits::gamma[index], o->GetGamma(), ntau*sizeof(slibreal_t));

    Orbits::classification[index] = o->GetClassification();

    Orbits::drift_shift[index] = o->GetDriftShift();
	Orbits::f[index] = part->GetF();
}

/**
 * Initialize this tool before using it (but after
 * configuration).
 */
void Orbits::Initialize() { }

/**
 * Enables or disables computation of the Jacobian's related
 * to the orbit.
 */
void Orbits::ToggleJacobianCalculation(bool compJ) {
    this->computeJacobian = compJ;
    pusher->ToggleJacobianCalculation(compJ);
}

