#ifndef _ANGULAR_DISTRIBUTION_H
#define _ANGULAR_DISTRIBUTION_H

#include <softlib/config.h>
#include "SOFT.h"
#include "Tools/Radiation/Models/Model.h"
#include "Tools/Radiation/Models/AngularDistribution/ADEmission.h"
#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"
#include "Tools/Radiation/Models/AngularDistribution/ADQuadrature2D.h"
#include "Tools/Radiation/RadiationException.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class AngularDistribution : public Model {
        private:
            ADEmission *emission;
            ADQuadrature2D *quadrature2d;
            slibreal_t power,
                *I, *Q, *U, *V;
            
            // General variable used to store prefactors
            slibreal_t prefactor;

            // Number of points on detector surface
            unsigned int nsamples;

            std::string emissionName = "N/A";

            // Default settings
            static const unsigned int NSAMPLES_DEFAULT;
        public:
            AngularDistribution(Radiation *rad);
            ~AngularDistribution();

            virtual void Configure(struct global_settings*, ConfigBlock*, ConfigBlock*) override;
            void ConfigureEmission(ConfigBlock*, const std::string&, struct global_settings*);
            ADSynchrotronEmission *ConfigureSynchrotronEmission(ConfigBlock*, struct global_settings*);
            virtual void InitializeOrbit(Orbit *__UNUSED__(o)) override {}
            virtual void InitializeTimestep(RadiationParticle*) override;

            virtual void HandleParticle(RadiationParticle*, const slibreal_t, const slibreal_t) override;

            virtual const std::string GetDescription() const override;

            virtual slibreal_t GetPower() override { return this->emission->GetTotalEmission(); }
            virtual slibreal_t GetPowerQ() override { return this->emission->GetTotalEmissionQ(); }
            virtual slibreal_t GetPowerU() override { return this->emission->GetTotalEmissionU(); }
            virtual slibreal_t GetPowerV() override { return this->emission->GetTotalEmissionV(); }
            virtual slibreal_t *GetSpectrum() override { return this->I; }
            virtual slibreal_t *GetStokesI() override { return this->I; }
            virtual slibreal_t *GetStokesQ() override { return this->Q; }
            virtual slibreal_t *GetStokesU() override { return this->U; }
            virtual slibreal_t *GetStokesV() override { return this->V; }
            virtual slibreal_t *GetWavelengths() override { return this->emission->GetWavelengths(); }
            virtual unsigned int GetNWavelengths() override { return this->emission->GetNWavelengths(); }
    };
}

#endif/*_ANGULAR_DISTRIBUTION_H*/
