#ifndef _CONE_H
#define _CONE_H

#include <string>
#include <softlib/config.h>
#include "SOFT.h"
#include "Tools/Radiation/RadiationException.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/Model.h"
#include "Tools/Radiation/Models/Cone/ConeEmission.h"
#include "Tools/Radiation/Models/Cone/Projection/ConeProjection.h"

namespace __Radiation {
    class Cone : public Model {
        private:
            ConeEmission *emission;
            ConeProjection *projection;

            bool edgeCheck;
            bool precomputeEmission;
            slibreal_t overlapFraction = 1;
            
            slibreal_t totEmission,
                totQ, totU, totV,
                *I, *Q, *U, *V;

            std::string emissionName = "N/A";
        public:
            Cone(Radiation *rad);
            ~Cone();

            virtual void Configure(struct global_settings*, ConfigBlock*, ConfigBlock*) override;
            void ConfigureEmission(const std::string&, ConfigBlock*);
            virtual void InitializeOrbit(Orbit *__UNUSED__(o)) override {}
            virtual void InitializeTimestep(RadiationParticle*) override;

            virtual void HandleParticle(RadiationParticle*, const slibreal_t, const slibreal_t) override;

            virtual const std::string GetDescription() const override;

            virtual slibreal_t GetPower() override { return this->totEmission; }
            virtual slibreal_t GetPowerQ() override { return this->totQ; }
            virtual slibreal_t GetPowerU() override { return this->totU; }
            virtual slibreal_t GetPowerV() override { return this->totV; }
            virtual slibreal_t *GetSpectrum() override { return this->I; }
            virtual slibreal_t *GetStokesI() override { return this->I; }
            virtual slibreal_t *GetStokesQ() override { return this->Q; }
            virtual slibreal_t *GetStokesU() override { return this->U; }
            virtual slibreal_t *GetStokesV() override { return this->V; }
            virtual slibreal_t *GetWavelengths() override { return this->emission->GetWavelengths(); }
            virtual unsigned int GetNWavelengths() override { return this->emission->GetNWavelengths(); }

            bool EdgeCheck(__Radiation::RadiationParticle*, const slibreal_t, const slibreal_t);
    };

    class ConeException : public RadiationException {
        public:
            template<typename ... Args>
            ConeException(const std::string &msg, Args&& ... args)
                : RadiationException(msg, std::forward<Args>(args) ...) {
                AddModule("Cone model");
            }
    };
}

#endif/*_CONE_H*/
