#ifndef _ORBIT_H
#define _ORBIT_H

class Orbit;

typedef enum {
	ORBIT_CLASS_UNKNOWN,
	ORBIT_CLASS_PASSING,
	ORBIT_CLASS_TRAPPED
} orbit_class_t;

#include <softlib/config.h>
#include <softlib/Integrator.h>
#include "PhaseSpace/Particle.h"
#include "Orbit/SOFTEquation.h"

class Orbit {
	protected:
		// Classification of orbit topology
		orbit_class_t orbitClass = ORBIT_CLASS_UNKNOWN;
		// Whether orbit was completed (false if particle collided with wall)
		bool completed = false;
        // True if magnetic field derivatives have been evaluated along this orbit
        bool hasBDerivatives = false;

		unsigned int ntau = 0;
        // Phase space indices
        unsigned int ir, ip1, ip2;
        slibreal_t m, q;       // Particle mass and charge

		/* Time-evolving quantities */
		slibreal_t
			*tau,		/* Time points (orbit time) */
			*x,			/* Particle/guiding-center position vector */
			*p,			/* Particle/guiding-center momentum vector */
			*ppar,		/* Parallel momentum of particle */
			*pperp,		/* Perpendicular momentum of particle */
			*Jdtdrho,	/* Trajectory coordinate Jacobian (times dtau * drho) */
            *_solution, /* Temporary storage of integrator solution */
            *_solution2;/* Temporary storage of secondary integrator solution */

        /* Related vectors */
        slibreal_t
            *Babs,          /* Magnetic field strength vector (1-dimensional) */
            *B,             /* Magnetic field vector (3-dimensional) */
            *bhat,          /* Magnetic field unit vector (3-dimensional) */
            *p2,            /* Momentum squared */
            *ppar2,         /* Parallel momentum squared */
            *pperp2,        /* Perpendicular momentum squared */
            *gamma;         /* Lorentz factor (or energy) */

        /* Vectors only stored when including drift orbits */
        slibreal_t
            *gradB,         /* Gradient of magnetic field strength */
            *curlB,         /* Curl of magnetic field vector */
            **jacobianB;    /* Magnetic field jacobian */
	public:
        Orbit(unsigned int, bool calcBDerivatives=false);
        ~Orbit();
        void CopyTo(Orbit*);
		Orbit *Create(slibreal_t, Integrator<6>*, Integrator<6>*, SOFTEquation*, Particle*, slibreal_t, orbit_class_t);

        unsigned int GetNTau() const { return this->ntau; }
        bool HasBDerivatives() const { return this->hasBDerivatives; }

        slibreal_t *GetTau() const { return this->tau; }
        slibreal_t GetTau(unsigned int ti) const { return this->tau[ti]; }
        slibreal_t *GetX() const { return this->x; }
        slibreal_t *GetX(unsigned int ti) const { return (this->x+3*ti); }
        slibreal_t *GetP() const { return this->p; }
        slibreal_t *GetP(unsigned int ti) const { return (this->p+3*ti); }
        slibreal_t *GetPpar() const { return this->ppar; }
        slibreal_t GetPpar(unsigned int ti) const { return this->ppar[ti]; }
        slibreal_t *GetPperp() const { return this->pperp; }
        slibreal_t GetPperp(unsigned int ti) const { return this->pperp[ti]; }
        slibreal_t *GetJdtdrho() const { return this->Jdtdrho; }
        slibreal_t GetJdtdrho(unsigned int ti) const { return this->Jdtdrho[ti]; }

        orbit_class_t GetClassification() const { return this->orbitClass; }
        void SetClassification(const orbit_class_t v) { this->orbitClass = v; }

        slibreal_t *GetInternalSolution() const { return this->_solution; }
        slibreal_t *GetInternalSolution(unsigned int ti) const { return (this->_solution+6*ti); }
        slibreal_t *GetInternalSolutionSecondary() const { return this->_solution2; }
        slibreal_t *GetInternalSolutionSecondary(unsigned int ti) const { return (this->_solution2+6*ti); }

        slibreal_t *GetBabs() const { return this->Babs; }
        slibreal_t GetBabs(unsigned int ti) const { return this->Babs[ti]; }
        slibreal_t *GetB() const { return this->B; }
        slibreal_t *GetB(unsigned int ti) const { return (this->B+3*ti); }
        slibreal_t *GetBhat() const { return this->bhat; }
        slibreal_t *GetBhat(unsigned int ti) const { return (this->bhat+3*ti); }
        slibreal_t *GetGradB() const { if (hasBDerivatives) return this->gradB; else return nullptr; }
        slibreal_t *GetGradB(unsigned int ti) const { if (hasBDerivatives) return this->gradB+3*ti; else return nullptr; }
        slibreal_t *GetCurlB() const { if (hasBDerivatives) return this->curlB; else return nullptr; }
        slibreal_t *GetCurlB(unsigned int ti) const { if (hasBDerivatives) return this->curlB+3*ti; else return nullptr; }
        slibreal_t **GetBJacobian() const { if (hasBDerivatives) return this->jacobianB; else return nullptr; }
        slibreal_t **GetBJacobian(unsigned int ti) const { if (hasBDerivatives) return this->jacobianB+3*ti; else return nullptr; }
        slibreal_t *GetP2() const { return this->p2; }
        slibreal_t GetP2(unsigned int ti) const { return this->p2[ti]; }
        slibreal_t *GetPpar2() const { return this->ppar2; }
        slibreal_t GetPpar2(unsigned int ti) const { return this->ppar2[ti]; }
        slibreal_t *GetPperp2() const { return this->pperp2; }
        slibreal_t GetPperp2(unsigned int ti) const { return this->pperp2[ti]; }
        slibreal_t *GetGamma() const { return this->gamma; }
        slibreal_t GetGamma(unsigned int ti) const { return this->gamma[ti]; }

        slibreal_t GetCharge() const { return this->q; }
        slibreal_t GetMass() const { return this->m; }
        unsigned int GetIndexR() const { return this->ir; }
        unsigned int GetIndexP1() const { return this->ip1; }
        unsigned int GetIndexP2() const { return this->ip2; }
};

class OrbitException : public SOFTException {
    public:
        template<typename ... Args>
        OrbitException(const std::string &msg, Args&& ... args)
            : SOFTException(msg, std::forward<Args>(args) ...) {
            AddModule("Orbit");
        }
};

#endif/*_ORBIT_H*/
