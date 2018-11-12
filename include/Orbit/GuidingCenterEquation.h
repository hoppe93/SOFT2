#ifndef _GUIDING_CENTER_EQUATION_H
#define _GUIDING_CENTER_EQUATION_H

#include <omp.h>
#include <softlib/IntegratorEquation.h>
#include <softlib/MagneticField/MagneticField2D.h>

#include "SOFT.h"
#include "Orbit/SOFTEquation.h"

class GuidingCenterEquation : public SOFTEquation {
	private:
		bool include_drifts=false;
		SOFT *soft;
	public:
        static const int
            COORD_X   =0,
            COORD_Y   =1,
            COORD_Z   =2,
            COORD_PPAR=3,
            COORD_MU  =4,
            COORD_ZETA=5;

        GuidingCenterEquation(MagneticField2D *mf, struct global_settings *gs)
            : SOFTEquation(mf, gs) {
            this->include_drifts = gs->include_drifts;
        }

        virtual orbit_class_t ClassifyOrbit(Integrator<6>*);
        orbit_class_t ClassifyOrbitPpar(const slibreal_t*, const unsigned int);

        Vector<6>& Evaluate(const slibreal_t T, const Vector<6>& zval, Vector<6>& dzdt) { return Evaluate(T,zval,dzdt,nullptr); }
		Vector<6>& Evaluate(const slibreal_t, const Vector<6>&, Vector<6>&, slibreal_t* gamma);
        slibreal_t GetPositionR(slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t);
        slibreal_t GetPositionZ(slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t);
        bool IncludesDrifts() { return this->include_drifts; }
        Vector<6>& InitializeParticle(Particle*, Vector<6>&);
        void ToOrbitQuantities(slibreal_t*, slibreal_t*, Orbit*, slibreal_t, orbit_class_t);
		void ToggleDrifts(bool);
};

#endif/*_GUIDING_CENTER_EQUATION_H*/
