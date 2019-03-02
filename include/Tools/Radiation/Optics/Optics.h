#ifndef _RADIATION_OPTICS_H
#define _RADIATION_OPTICS_H

#include <complex>
#include <softlib/Vector.h>

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
                std::complex<slibreal_t> *Ex, *Ey;
                Vector<3> xhat, yhat, zhat;
            };

            // Pointer to simulation detector
            Detector *detector;
        private:
        public:
            Optics(Detector *det) { this->detector = det; }

            virtual void ApplyOptics(
                const struct Efield&,
                slibreal_t *I, slibreal_t *Q=nullptr,
                slibreal_t *U=nullptr, slibreal_t *V=nullptr
            ) = 0;
    };
}

#endif/*_RADIATION_OPTICS_H*/
