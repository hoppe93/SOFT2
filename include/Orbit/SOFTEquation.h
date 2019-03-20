#ifndef _SOFT_EQUATION_H
#define _SOFT_EQUATION_H

#include <softlib/IntegratorEquation.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/Vector.h>

class SOFTEquation;

#include "Orbit/Orbit.h"
#include "PhaseSpace/Particle.h"
#include "SOFT.h"

class SOFTEquation : public IntegratorEquation<6> {
    protected:
        MagneticField2D *magfield;
        struct global_settings *globset;
        Particle *particle;
    public:
        SOFTEquation(MagneticField2D *mf, struct global_settings *gs) {
            this->magfield = mf;
            this->globset = gs;
        }

        void CalculateJacobians(slibreal_t*, slibreal_t*, Orbit*, slibreal_t, bool);
        virtual orbit_class_t ClassifyOrbit(Integrator<6>*) = 0;

        virtual slibreal_t GetPositionR(
            slibreal_t, slibreal_t, slibreal_t,
            slibreal_t, slibreal_t, slibreal_t
        ) = 0;
        slibreal_t GetPositionR(const slibreal_t*);
        slibreal_t GetPositionR(const Vector<6>&);
        virtual slibreal_t GetPositionZ(
            slibreal_t, slibreal_t, slibreal_t,
            slibreal_t, slibreal_t, slibreal_t
        ) = 0;
        slibreal_t GetPositionZ(const slibreal_t*);
        slibreal_t GetPositionZ(const Vector<6>&);

        virtual Vector<6>& InitializeParticle(Particle*, Vector<6>&) = 0;
        virtual void ToOrbitQuantities(slibreal_t*, slibreal_t*, Orbit*, slibreal_t, orbit_class_t, bool) = 0;
};

#endif/*_SOFT_EQUATION*/
