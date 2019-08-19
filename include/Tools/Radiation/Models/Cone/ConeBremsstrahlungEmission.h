#ifndef _CONE_BREMSSTRAHLUNG_EMISSION_H
#define _CONE_BREMSSTRAHLUNG_EMISSION_H

#include <softlib/constants.h>
#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Cone/ConeEmission.h"
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

slibreal_t dilog_func(const slibreal_t);

namespace __Radiation {
    class ConeBremsstrahlungEmission : public ConeEmission {
        private:
#define _c  LIGHTSPEED
#define _c2 (_c*_c)
#define _e  ELEMENTARY_CHARGE
#define _e2 (_e*_e)
#define _m  ELECTRON_MASS
#define _m2 (_m*_m)
            static constexpr slibreal_t alpha =
                _e*_e / (4.0*M_PI*EPS0*_c*HBAR);
            static constexpr slibreal_t r02Alpha = 
                _e2*_e2*alpha / (16.0*M_PI*M_PI*EPS0*EPS0*_m2*_c2*_c2);

            unsigned int nspecies;
            slibreal_t *Z;
            slibreal_t *density;
        public:
            ConeBremsstrahlungEmission(Detector *det, MagneticField2D *mf, unsigned int nspecies, slibreal_t *Z, slibreal_t *density);
            ~ConeBremsstrahlungEmission();

            void HandleParticle(RadiationParticle*, bool);

            void CalculatePolarization(RadiationParticle*);
            void CalculateSpectrum(RadiationParticle*);
            void CalculateTotalEmission(RadiationParticle*);

            void IntegrateSpectrum();
    };

    class ConeBremsstrahlungException : public ConeException {
        public:
            template<typename ... Args>
            ConeBremsstrahlungException(const std::string &msg, Args&& ... args)
                : ConeException(msg, std::forward<Args>(args) ...) {
                AddModule("ConeBremsstrahlung");
            }
    };
}

#endif/*_CONE_BREMSSTRAHLUNG_EMISSION_H*/
