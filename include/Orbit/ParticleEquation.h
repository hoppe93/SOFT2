#ifndef _PARTICLE_EQUATION_H
#define _PARTICLE_EQUATION_H

#include <softlib/IntegratorEquation.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "PhaseSpace/Particle.h"
#include "Orbit/SOFTEquation.h"

class ParticleEquation : public SOFTEquation {
	public:
		ParticleEquation(MagneticField2D *mf, struct global_settings *gs)
            : SOFTEquation(mf, gs) {}

        virtual orbit_class_t ClassifyOrbit(Integrator<6>*);
        orbit_class_t ClassifyOrbitPpar(const slibreal_t*, const unsigned int);

        virtual orbit_type_t GetOrbitType() const { return ORBIT_TYPE_PARTICLE; }

		Vector<6>& Evaluate(const slibreal_t, const Vector<6>&, Vector<6>&);
        slibreal_t GetPositionR(slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t);
        slibreal_t GetPositionZ(slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t);
		Vector<6>& InitializeParticle(Particle*, Vector<6>&);
        void ToOrbitQuantities(slibreal_t*, slibreal_t*, slibreal_t*, Orbit*, slibreal_t, orbit_class_t, bool);
};

#endif/*_PARTICLE_EQUATION_H*/
