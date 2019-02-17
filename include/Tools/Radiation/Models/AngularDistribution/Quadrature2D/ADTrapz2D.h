#ifndef _ANGULAR_DISTRIBUTION_TRAPZ_2D_H
#define _ANGULAR_DISTRIBUTION_TRAPZ_2D_H

#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/Models/AngularDistribution/ADEmission.h"
#include "Tools/Radiation/Models/AngularDistribution/ADQuadrature2D.h"

namespace __Radiation {
    class ADTrapz2D : public ADQuadrature2D {
        private:
            unsigned int nsamples;
        public:
            ADTrapz2D(ADEmission*, Detector*, unsigned int);

            slibreal_t Integrate(
                RadiationParticle*, bool,
                slibreal_t*, slibreal_t*,
                slibreal_t*, slibreal_t*
            );

            template<bool, bool>
            slibreal_t _OuterIntegral(
                RadiationParticle*,
                slibreal_t*,
                slibreal_t *Q=nullptr, slibreal_t *U=nullptr, slibreal_t *V=nullptr
            );

            template<bool, bool>
            slibreal_t _InnerIntegral(
                slibreal_t, slibreal_t, slibreal_t,
                RadiationParticle*
            );
    };
}

#endif/*_ANGULAR_DISTRIBUTION_TRAPZ_2D_H*/
