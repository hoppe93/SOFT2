#ifndef _CONE_BREMSSTRAHLUNG_SCREENED_EMISSION_H
#define _CONE_BREMSSTRAHLUNG_SCREENED_EMISSION_H

#include <softlib/constants.h>
#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Cone/ConeEmission.h"
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

#include <gsl/gsl_integration.h>

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
            slibreal_t r02;

            double qagsEpsAbs=0.0,
                   qagsEpsRel=1e-3;
            std::size_t qagsLimit = 100;
            gsl_integration_workspace *qagsWS = nullptr;

       public:
            ConeBremsstrahlungScreenedEmission(Detector *det, MagneticField2D *mf, unsigned int nspecies, slibreal_t *Z, slibreal_t *Z0, slibreal_t *density);
            ~ConeBremsstrahlungScreenedEmission();

            void HandleParticle(RadiationParticle*, bool);

            void CalculatePolarization(RadiationParticle*);

            struct func_params {slibreal_t Z; slibreal_t Nej; slibreal_t a_bar; slibreal_t q0; slibreal_t lnq0;}; //Data-type for function parameters 
            void CalculateSpectrum(RadiationParticle*);
            slibreal_t FirstSpectrumIntegral(slibreal_t, slibreal_t, slibreal_t, slibreal_t); 
            slibreal_t SecondSpectrumIntegral(slibreal_t, slibreal_t, slibreal_t, slibreal_t);
            static double First_Integrand(slibreal_t, void*);
            static double Second_Integrand(slibreal_t, void*);
            static double Second_Integrand_app_sq0(slibreal_t, void*); 
            static double Second_Integrand_app_lq0(slibreal_t, void*); 
            static slibreal_t CalculateFormFactor(slibreal_t, slibreal_t, slibreal_t);
            
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
