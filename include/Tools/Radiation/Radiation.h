#ifndef _RADIATION_H
#define _RADIATION_H

namespace __Radiation {
    class Radiation;
}

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/DistributionFunction/DistributionFunction.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/Vector.h>
#include "Orbit/Orbit.h"
#include "PhaseSpace/Particle.h"
#include "PhaseSpace/ParticleGenerator.h"
#include "SOFT.h"
#include "Tools/Tool.h"
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/Models/Model.h"
#include "Tools/Radiation/Output/RadiationOutput.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    struct radiation_modelspec {
        std::string name;
        Model *(*init)(struct global_settings*, ConfigBlock*, ConfigBlock*, Radiation*);
    };
    struct radiation_outputspec {
        std::string name;
        RadiationOutput *(*init)(ConfigBlock*, ConfigBlock*, Detector*, MagneticField2D*, ParticleGenerator*);
    };

    // Stokes parameter identifiers
    typedef enum {
        STOKES_PARAM_I,
        STOKES_PARAM_Q,
        STOKES_PARAM_U,
        STOKES_PARAM_V
    } stokes_param_t;

    typedef enum {
        WALL_OPACITY_OPAQUE,            // Full wall checks
        WALL_OPACITY_SEMI_TRANSPARENT,  // No outer wall
        WALL_OPACITY_TRANSPARENT        // No wall checks
    } wall_opacity_t;

    class Radiation : public Tool {
        public:
            enum RadiationQuadrature {
                QUADRATURE_FINDSOV,
                QUADRATURE_TRAPEZOIDAL
            };

            static int CONFBLOCK_T_DETECTOR,
                       CONFBLOCK_T_DETECTOR_OPTICS,
                       CONFBLOCK_T_MODEL,
                       CONFBLOCK_T_OUTPUT;
        private:
            // Toroidal resolution
            unsigned int ntoroidal;
            // Toroidal step increments
            //slibreal_t dphi, cosdphi, sindphi;
            slibreal_t dphi, *cosphi, *sinphi;

            enum RadiationQuadrature quadrature = QUADRATURE_FINDSOV;
            char *torflags;
            slibreal_t torthreshold = 1e-9;

            bool ignoreTrapped = false;
            bool measuresPolarization = false;
            bool shiftLarmorRadius = false;
            wall_opacity_t wall_opacity;
        public:
            __Radiation::Detector *detector;
            __Radiation::Model *model;

            MagneticField2D *magfield;
            ParticleGenerator *pgen;

            unsigned int noutput = 0;
            __Radiation::RadiationOutput **output = nullptr;

            Radiation() : Tool("Radiation") {};
            Radiation(MagneticField2D*, ParticleGenerator*, ParticlePusher*);
            ~Radiation();

            virtual void Configure(struct global_settings*, ConfigBlock*, ConfigBlock*) override;
            virtual void Finish() override;
            virtual void Handle(Orbit*, Particle*) override;
            virtual void Initialize() override;
            virtual void Output() override;
            void RegisterOutput(RadiationParticle*);
            Model *SetupRadiationModel(struct global_settings*, ConfigBlock*, ConfigBlock*);
            RadiationOutput *SetupRadiationOutput(ConfigBlock*, ConfigBlock*);
            virtual void Welcome(const std::string &prefix="  ") override;

            void HandleTimeIntegral(
                Orbit*, Particle*, void (Radiation::*)(
                    RadiationParticle&, orbit_type_t, slibreal_t,
                    slibreal_t, slibreal_t, slibreal_t, slibreal_t
                )
            );
            void EvaluateToroidalTrapz(
                RadiationParticle&, orbit_type_t,
                slibreal_t, slibreal_t, slibreal_t,
                slibreal_t, slibreal_t
            );
            void EvaluateToroidalTrapzImproved(
                RadiationParticle&, orbit_type_t,
                slibreal_t, slibreal_t, slibreal_t,
                slibreal_t, slibreal_t
            );
            // Regular trapezoidal rule
            void HandleTrapz(Orbit*, Particle*);
            // Improved trapezoidal rule
            void HandleTrapzImproved(Orbit*, Particle*);
            unsigned int IntegrateToroidalImproved(RadiationParticle&, orbit_type_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t, int, int, slibreal_t&);
            void LocateSurfaceOfVisibility(RadiationParticle*, unsigned int*, unsigned int*);
            unsigned int LocatePointOfVisibility(RadiationParticle*);
            void ResetToroidalFlags();

            void SetDetector(Detector *det) { this->detector = det; }
            void SetModel(Model *mdl) { this->model = mdl; }

            bool IsWithinFieldOfView(slibreal_t, slibreal_t, slibreal_t, Vector<3>&);
            bool MeasuresPolarization() { return measuresPolarization; }

            void SetShiftedLarmorRadius(bool v) { this->shiftLarmorRadius = v; }
            void ShiftLarmorRadius(slibreal_t&, slibreal_t&, slibreal_t&, const slibreal_t, const slibreal_t, const slibreal_t, RadiationParticle*);

            static void PrepareConfiguration(Configuration*);
    };
}

#endif/*_RADIATION_H*/
