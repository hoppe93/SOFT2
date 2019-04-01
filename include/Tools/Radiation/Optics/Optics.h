#ifndef _RADIATION_OPTICS_H
#define _RADIATION_OPTICS_H

#include <softlib/config.h>
#include <softlib/Vector.h>
#include "Tools/Radiation/Detector.h"

namespace __Radiation {
    class Optics {
        public:
            /**
             * This Efield structure is intended to contain
             * multiple frequency components of the electric
             * field emitted in the 'zhat' direction. The 'xhat'
             * and 'yhat' vectors indicate the two possible
             * (orthonormal) polarization directions of the
             * electric field.
             */
            struct Efield {
                unsigned int nE;
                slibreal_t *Ex2, *Ey2, *ExEy;
                Vector<3> xhat, yhat, zhat;
            };
        protected:
            // Pointer to simulation detector
            Detector *detector;
        public:
            Optics(Detector *det) { this->detector = det; }

            virtual void ApplyOptics(
                const struct Optics::Efield&,
                slibreal_t *I, slibreal_t *Q=nullptr,
                slibreal_t *U=nullptr, slibreal_t *V=nullptr
            ) = 0;
    };
}

#endif/*_RADIATION_OPTICS_H*/
