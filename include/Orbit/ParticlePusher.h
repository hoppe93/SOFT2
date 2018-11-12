#ifndef _PARTICLE_PUSHER_H
#define _PARTICLE_PUSHER_H

#include <string>
#include <cmath>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/Integrator.h>
#include <softlib/IntegratorEquation.h>
#include <softlib/MagneticField/MagneticField2D.h>

class ParticlePusher;

#include "Orbit/Orbit.h"
#include "SOFT.h"
#include "SOFTEquation.h"

enum orbittimeunit {
    ORBITTIMEUNIT_POLOIDAL,
    ORBITTIMEUNIT_SECONDS
};

class ParticlePusher {
	private:
        bool calculateJacobianOrbit = true;
		ConfigBlock settings;
		Integrator<6> *integrator1, *integrator2;
        MagneticField2D *magfield;
        struct global_settings* globset;
        SOFTEquation *equation;

        slibreal_t integrator_tol;
        slibreal_t maxtime;
        enum orbittimeunit timeunit;

        // The value by which a particle is 'nudged' in the
        // radial direction in order to solve for a second
        // orbit which can be used to calculate the spatial
        // Jacobian determinant.
        slibreal_t nudge_value;

        Orbit *retorbit;
        slibreal_t *solution;           // Temporary storage of 6D solution

        // Previous and initial R and Z position of particle
        // (for calculating duration of a poloidal orbit)
        slibreal_t rprev, zprev,
                   rinit, zinit;
        unsigned int zpass,             // Number of times particle passed z = z0
                     restflag;          // Number of timesteps during which the particle
                                        // hasn't moved vertically.

        uint32_t ntimesteps;            // Number of time steps per orbit

        Orbit *CreateOrbit(Orbit*, slibreal_t);

	public:
        ParticlePusher(MagneticField2D*, struct global_settings *globset);
		ParticlePusher(MagneticField2D*, struct global_settings *globset, ConfigBlock*, ConfigBlock*);

        static const std::string equation_defaults;

        void EvaluateSecondaryOrbit(Particle*, enum Particle::nudge_direction);
		void InitDefaults();
        void InitEquation(const std::string&, ConfigBlock&);
        void InitGeneralIntegrator(ConfigBlock&, IntegratorEquation<6>*);

        slibreal_t FindPoloidalTime(SOFTEquation*, Integrator<6>*);
        uint32_t GetNTimeSteps() { return ntimesteps; }
        slibreal_t GetPoloidalTime(slibreal_t);
        slibreal_t GetPositionR(SOFTEquation*, Integrator<6>*);
        slibreal_t GetPositionZ(SOFTEquation*, Integrator<6>*);
        bool MaxTimeReached(slibreal_t, slibreal_t);
        void ResetPoloidalTime();

        /* Mainly for internal use */
        Integrator<6> *GetIntegrator1() { return integrator1; }
        Orbit *GetOrbit() { return retorbit; }

        Orbit *Push(Particle*);
        void RunIntegrator(SOFTEquation*, Integrator<6>*, Particle*);

        void ToggleJacobianCalculation(bool);
};

class ParticlePusherException : public SOFTException {
    public:
        template<typename ... Args>
        ParticlePusherException(const std::string &msg, Args&& ... args)
            : SOFTException(msg, std::forward<Args>(args) ...) {
            AddModule("Particle pusher");
        }
};

#endif/*_PARTICLE_PUSHER_H*/
