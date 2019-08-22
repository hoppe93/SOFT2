#ifndef _ISOTROPIC_H
#define _ISOTROPIC_H

#include <string>
#include <softlib/config.h>
#include "SOFT.h"
#include "Tools/Radiation/RadiationException.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/Model.h"

namespace __Radiation {
    class Isotropic : public Model {
        private:
            slibreal_t value, power;
        public:
            Isotropic(Radiation *rad) : Model(rad) {nonzero=true;}

            virtual void Configure(struct global_settings*, ConfigBlock*, ConfigBlock*) override;
            virtual void InitializeOrbit(Orbit *__UNUSED__(o)) override {};
            virtual void InitializeTimestep(RadiationParticle *__UNUSED__(rp)) override {};

            virtual void HandleParticle(RadiationParticle*, orbit_type_t, const slibreal_t, const slibreal_t) override;

            virtual const std::string GetDescription() const override;

            static slibreal_t LambdaIso(slibreal_t, slibreal_t, slibreal_t, slibreal_t, slibreal_t);

            void SetRadiationValue(slibreal_t);

            virtual slibreal_t GetPower() override { return this->power; }
            virtual slibreal_t GetPowerQ() override { return 0; }
            virtual slibreal_t GetPowerU() override { return 0; }
            virtual slibreal_t GetPowerV() override { return 0; }
            virtual slibreal_t *GetSpectrum() override { return nullptr; }
            virtual slibreal_t *GetStokesI() override { return nullptr; }
            virtual slibreal_t *GetStokesQ() override { return nullptr; }
            virtual slibreal_t *GetStokesU() override { return nullptr; }
            virtual slibreal_t *GetStokesV() override { return nullptr; }
            virtual slibreal_t *GetWavelengths() override { return nullptr; }
            virtual unsigned int GetNWavelengths() override { return 0; }
    };

    class IsotropicException : public RadiationException {
        public:
            template<typename ... Args>
            IsotropicException(const std::string &msg, Args&& ... args)
                : RadiationException(msg, std::forward<Args>(args) ...) {
                AddModule("Isotropic");
            }
    };
}

#endif/*_ISOTROPIC_H*/
