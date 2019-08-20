#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <softlib/config.h>
#include <softlib/Vector.h>
#include "constants.h"
#include "SOFTException.h"

class Particle {
    public:
        enum nudge_direction { NUDGE_INWARDS, NUDGE_OUTWARDS };
		enum coordinate {
			COORDINATE_GAMMA,
			COORDINATE_P,
			COORDINATE_PPAR,
			COORDINATE_PPERP,
			COORDINATE_THETAP,
			COORDINATE_ITHETAP,
			COORDINATE_XI
		};
		enum position {
			POSITION_GUIDINGCENTER,
			POSITION_PARTICLE
		};

	private:
		slibreal_t rho,					// Particle radius (trajectory coordinate)
                   z0;                  // Initial vertical position
		slibreal_t drift_shift;			// Drift orbit shift (in 1st order GC theory)
		enum position position_type;	// Which position is specified (particle or GC).
										// Values according to 'Particle::POSITION_XXX' below.

		slibreal_t pmag;				// Magnitude of p
		slibreal_t ppar, pperp, zeta;	// Particle momentum in cylindrical coordinates
		slibreal_t thetap, xi;			// Pitch angle and cosine of pitch angle

        slibreal_t f;                   // Value of distribution function

        // Phase-space indices corresponding
        // to this particle
        unsigned int ir, ip1, ip2, izeta;

		// Useful quantities
		slibreal_t gamma;				// Particle Lorentz factor
	
		enum coordinate momentum1,      // Types of momentum coordinates (values
                        momentum2;	    // according to 'Particle::COORDINATE_XXX' below)
										

		slibreal_t drho,				// Differential elements for spatial and momentum-space coordinates
				dparam1, dparam2, dzeta;		

		slibreal_t m=ELECTRON_MASS,     // Particle mass
                   q=-ELEMENTARY_CHARGE;// Particle charge

		bool momentum_initialized=false,// True when the two momentum coordinates have been initialized
			 position_initialized=false;// True when rho has been set

		/* Methods for initializing all momentum coordinates
		 * (except cartesian) from two coordinates. */
		void Init_GammaPPar(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		void Init_GammaThetap(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		void Init_GammaXi(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		void Init_PPPar(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		void Init_PThetap(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		void Init_PXi(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		void Init_PParPPerp(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		void Init_PParThetap(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		void Init_PParXi(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		void Init_PPerpThetap(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		void Init_PPerpXi(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);

        /* Methods for quickly converting a pair of momentum
         * coordinates to a ppar/pperp pair. */
		static void ToPP_GammaPPar(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		static void ToPP_GammaThetap(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		static void ToPP_GammaXi(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		static void ToPP_PPPar(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		static void ToPP_PThetap(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		static void ToPP_PXi(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		static void ToPP_PParPPerp(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		static void ToPP_PParThetap(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		static void ToPP_PParXi(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		static void ToPP_PPerpThetap(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
		static void ToPP_PPerpXi(const slibreal_t, const slibreal_t, slibreal_t*, slibreal_t*);
	public:
		Particle();
		void InitializeMomentum(const enum coordinate, const enum coordinate, const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t);
		void InitializePosition(const enum position, const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t);
        void Nudge(const slibreal_t, enum nudge_direction nd=NUDGE_OUTWARDS);

		static const char* GetCoordinateName(const int);

        Vector<3> Get3Momentum(Vector<3>&);
        slibreal_t GetCharge()          const {return this->q;}
        slibreal_t GetMass()            const {return this->m;}
        enum position GetPositionType() const {return this->position_type;}

		/* Get coordinates */
		slibreal_t GetJMomentum1()  const { return this->dparam1; }
		slibreal_t GetJMomentum2()  const { return this->dparam2; }
        slibreal_t GetJZeta()       const { return this->dzeta; }

        slibreal_t GetActualRho()   const { return (this->rho+this->drift_shift); }
		slibreal_t GetDriftShift()  const { return this->drift_shift; }
        slibreal_t GetF()           const { return this->f; }
		slibreal_t GetRho()         const { return this->rho; }
        slibreal_t GetDRho()        const { return this->drho; }
        slibreal_t GetZ0()          const { return this->z0; }

		slibreal_t GetGamma()       const { return this->gamma; }
		slibreal_t GetMomentum()    const { return this->pmag; }
		slibreal_t GetPpar()        const { return this->ppar; }
		slibreal_t GetPperp()       const { return this->pperp; }
		slibreal_t GetThetap()      const { return this->thetap; }
		slibreal_t GetXi()          const { return this->xi; }

        slibreal_t GetZeta()        const {return this->zeta;}

        unsigned int GetIndexR()    const { return this->ir; }
        unsigned int GetIndexP1()   const { return this->ip1; }
        unsigned int GetIndexP2()   const { return this->ip2; }
        unsigned int GetIndexZeta() const { return this->izeta; }

		/* Set coordinates */
		void SetDRho(const slibreal_t v)        { this->drho = v; }
		void SetDMomentum1(const slibreal_t v)  { this->dparam1 = v; }
		void SetDMomentum2(const slibreal_t v)  { this->dparam2 = v; }
        void SetDZeta(const slibreal_t v)       { this->dzeta = v; }

		void SetDriftShift(const slibreal_t v)  { this->drift_shift = v; }
		void SetRho(const slibreal_t v)         { this->rho = v; }

        void SetF(const slibreal_t v)           { this->f = v; }
		void SetGamma(const slibreal_t v)       { this->gamma = v; }
		void SetMomentum(const slibreal_t v)    { this->pmag = v; }
		void SetPpar(const slibreal_t v)        { this->ppar = v; }
		void SetPperp(const slibreal_t v)       { this->pperp = v; }
		void SetThetap(const slibreal_t v)      { this->thetap = v; }
		void SetXi(const slibreal_t v)          { this->xi = v; }

		/* Set properties */
		void SetCharge(const slibreal_t v)      { this->q = v; }
		void SetMass(const slibreal_t v)        { this->m = v; }
		void SetPositionType(enum position v) { this->position_type = v; }

        void SetIndexR(unsigned int i)          { this->ir = i; }
        void SetIndex1(unsigned int i)          { this->ip1 = i; }
        void SetIndex2(unsigned int i)          { this->ip2 = i; }
        void SetIndexZeta(unsigned int i)       { this->izeta = i; }
        void SetIndices(unsigned int ir, unsigned int i1, unsigned int i2, unsigned int izeta) {
            SetIndexR(ir);
            SetIndex1(i1);
            SetIndex2(i2);
            SetIndexZeta(izeta);
        }

        static void VerifyCoordinateSpecification(slibreal_t, slibreal_t, unsigned int, int);
        static void ToPP(const slibreal_t, const slibreal_t, const int, const int, slibreal_t*, slibreal_t*);
};

class ParticleException : public SOFTException {
    public:
        template<typename ... Args>
        ParticleException(const std::string &msg, Args&& ... args)
            : SOFTException(msg, std::forward<Args>(args) ...) {
            AddModule("Particle");
        }
};

#endif/*_PARTICLE_H*/
