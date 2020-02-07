#ifndef _ORBITS_H
#define _ORBITS_H

#include <string>
#include <softlib/Configuration.h>
#include "Orbit/Orbit.h"
#include "PhaseSpace/Particle.h"
#include "Tools/OutputModule.h"
#include "Tools/Tool.h"

class Orbits : public Tool, public OutputModule {
    private:
        unsigned int norbits;
        unsigned int nr, np1, np2;
        unsigned int ntau;
        bool computeJacobian=false;
        size_t allocatedBytes;

        std::string output;

		static slibreal_t
			**tau,		/* Time points (orbit time) */
			**x,		/* Particle/guiding-center position vector */
			**p,		/* Particle/guiding-center momentum vector */
			**ppar,		/* Parallel momentum of particle */
			**pperp,	/* Perpendicular momentum of particle */
			**Jdtdrho,	/* Trajectory coordinate Jacobian (times dtau * drho) */
            **solution; /* Temporary storage of integrator solution */

        /* Related vectors */
        static slibreal_t
            **Babs,      /* Magnetic field strength vector (1-dimensional) */
            **B,         /* Magnetic field vector (3-dimensional) */
            **bhat,      /* Magnetic field unit vector (3-dimensional) */
            **p2,        /* Momentum squared */
            **ppar2,     /* Parallel momentum squared */
            **pperp2,    /* Perpendicular momentum squared */
            **gamma;     /* Lorentz factor (or energy) */

        static orbit_class_t *classification;   /* Orbit classifications */

        static const std::string DEFAULT_QUANTITIES[];
        static const unsigned int NDEFAULT_QUANTITIES;

        ParticlePusher *pusher; // Reference to pusher object
    public:
        Orbits(MagneticField2D*, ParticleGenerator*, ParticlePusher*);

        slibreal_t **Allocate(const unsigned int, const unsigned int, const unsigned int);
        orbit_class_t *Allocate_class(const unsigned int);
        virtual void Configure(struct global_settings*, ConfigBlock*, ConfigBlock*) override;
        virtual void Finish() override {}
        virtual void Handle(Orbit*, Particle*) override;
        virtual void Initialize() override;
        virtual void Output() override;
        void ToggleJacobianCalculation(bool);
        virtual void Welcome(const std::string &prefix="  ") override;

        static void PrepareConfiguration(Configuration*);
};

class OrbitsException : public ToolException {
    public:
        template<typename ... Args>
        OrbitsException(const std::string &msg, Args&& ... args)
            : ToolException(msg, std::forward<Args>(args) ...) {
            AddModule("Orbits");
        }
};

#endif/*_ORBITS_H*/
