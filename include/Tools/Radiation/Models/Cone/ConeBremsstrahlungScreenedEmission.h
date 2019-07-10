#ifndef _CONE_BREMSSTRAHLUNG_SCREENED_EMISSION_H
#define _CONE_BREMSSTRAHLUNG_SCREENED_EMISSION_H

#include <softlib/constants.h>
#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Cone/ConeEmission.h"
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class ConeBremsstrahlungScreenedEmission : public ConeEmission {
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
                _e2*alpha / (16.0*M_PI*M_PI*EPS0*EPS0*_m2*_c2);
            unsigned int nspecies;
            slibreal_t *Z;
            slibreal_t *Z0;
            slibreal_t *density;
            const slibreal_t r02 = ELECTRON_CHARGE*ELECTRON_CHARGE*ELECTRON_CHARGE*ELECTRON_CHARGE/(16*M_PI*M_PI*EPS0*EPS0*LIGHTSPEED*LIGHTSPEED*LIGHTSPEED*LIGHTSPEED*ELECTRON_MASS*ELECTRON_MASS);
            slibreal_t Z2;

       public:
            ConeBremsstrahlungScreenedEmission(Detector *det, MagneticField2D *mf, unsigned int nspecies, slibreal_t *Z, slibreal_t *Z0, slibreal_t *density)
                : ConeEmission(det, mf), nspecies(nspecies), Z(Z), Z0(Z0), density(density) {
            };
            ~ConeBremsstrahlungScreenedEmission() {}//delete []this->Z; delete []this->Z0; delete []this->density;} //Having the destructor with the delete statements currently makes the test not work as it tries to delete pointers. Maybe they should not be there at all?

            void HandleParticle(RadiationParticle*, bool);

            void CalculatePolarization(RadiationParticle*);
            void CalculateSpectrum(RadiationParticle*);
            void CalculateTotalEmission(); 
            slibreal_t Calculate4BS(slibreal_t);

            void IntegrateSpectrum();
    };

    class ConeBremsstrahlungScreenedException : public ConeException {
        public:
            template<typename ... Args>
            ConeBremsstrahlungScreenedException(const std::string &msg, Args&& ... args)
                : ConeException(msg, std::forward<Args>(args) ...) {
                AddModule("ConeBremsstrahlungScreened");
            }
    };
}

#endif/*_CONE_BREMSSTRAHLUNG_SCREENED_EMISSION_H*/
