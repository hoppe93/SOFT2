#ifndef _ANGULAR_DISTRIBUTION_EVAL2D_H
#define _ANGULAR_DISTRIBUTION_EVAL2D_H

#include "Tools/Radiation/Models/AngularDistribution/ADQuadrature2D.h"

namespace __Radiation {
    class ADEval2D : public ADQuadrature2D {
        private:
        public:
            ADEval2D(ADEmission *a, Detector *d) : ADQuadrature2D(a, d) {}

            slibreal_t Integrate(
                RadiationParticle*, bool,
                slibreal_t*, slibreal_t*, slibreal_t*, slibreal_t*
            );
    };
}

#endif/*_ANGULAR_DISTRIBUTION_EVAL2D_H*/
