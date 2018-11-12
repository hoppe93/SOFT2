#ifndef _ANGULAR_DISTRIBUTION_QUADRATURE_2D_H
#define _ANGULAR_DISTRIBUTION_QUADRATURE_2D_H

#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/Models/AngularDistribution/ADEmission.h"

namespace __Radiation {
    class ADQuadrature2D {
        protected:
            ADEmission *emission;
            Detector *detector;

            unsigned int nwavelengths;
            slibreal_t *I, *Q, *U, *V;
        public:
            ADQuadrature2D(ADEmission*, Detector*);
            virtual ~ADQuadrature2D();

            struct angles {
                slibreal_t
                    sinMu,  cosMu,
                    sin2Mu, cos2Mu,
                    rcp2;
                Vector<3> n;
            };

            Detector *GetDetector() const { return detector; }
            ADEmission *GetEmissionModel() const { return emission; }

            // Implemented members
            void EvaluateAngles(slibreal_t, slibreal_t, RadiationParticle*, struct angles&);

            /**
             * Multiply the elements of I, Q, U and V by the given
             * scalar. The result is
             *
             *   I[i] = a * I[i]
             *   etc.
             *
             * If the final four elements are specified, this function
             * instead does the following:
             *
             *   Id[i] = a * Is[i]
             *
             * a:              Scalar to multiply the arrays by.
             * n:              Number of elements in the arrays.
             * I, Q, U, V:     Arrays to multiply.
             * Is, Qs, Us, Vs: If provided, I-V are used as destinations.
             *
             * withPolarization: If false, multiplies only I. Otherwise,
             *                   all four arrays are multiplied by a.
             */
            template<bool withPolarization>
            void MultiplySpectra(
                slibreal_t a, unsigned int n,
                slibreal_t *I, slibreal_t *Q, slibreal_t *U, slibreal_t *V
            ) {
                unsigned int i;
                for (i = 0; i < n; i++) {
                    I[i] *= a;
                    if (withPolarization) {
                        Q[i] *= a;
                        U[i] *= a;
                        V[i] *= a;
                    }
                }
            }
            template<bool withPolarization>
            void MultiplySpectra(
                slibreal_t a, unsigned int n,
                slibreal_t *I, slibreal_t *Q, slibreal_t *U, slibreal_t *V,
                slibreal_t *Is, slibreal_t *Qs, slibreal_t *Us, slibreal_t *Vs
            ) {
                unsigned int i;
                for (i = 0; i < n; i++) {
                    I[i] = a * Is[i];
                    if (withPolarization) {
                        Q[i] = a * Qs[i];
                        U[i] = a * Us[i];
                        V[i] = a * Vs[i];
                    }
                }
            }

            /**
             * Sum the elements of Is, Qs, Us and Vs, (s for source)
             * into the destination arrays Id, Qd, Ud, Vd, while multiplying
             * the sources with the scalar a. The result is
             *
             *   Id[i] = Id[i] + a * Is[i];
             *   etc.
             * 
             * a:              Scalar to multiply the source elements with.
             * n:              Number of elements in each of the arrays.
             * Id, Qd, Ud, Vd: Destination arrays.
             * Is, Qs, Us, Vs: Source arrays.
             *
             * withPolarization: If false, sums only I. Otherwise, sums
             *                   all Stokes parameters.
             */
            template<bool withPolarization>
            void SumSpectra(
                slibreal_t a, unsigned int n,
                slibreal_t *Id, slibreal_t *Qd, slibreal_t *Ud, slibreal_t *Vd,
                slibreal_t *Is, slibreal_t *Qs, slibreal_t *Us, slibreal_t *Vs
            ) {
                unsigned int i;
                for (i = 0; i < n; i++) {
                    Id[i] += a * Is[i];
                    if (withPolarization) {
                        Qd[i] += a * Qs[i];
                        Ud[i] += a * Us[i];
                        Vd[i] += a * Vs[i];
                    }
                }
            }

            /**
             * Set the elements of the given arrays to zero.
             * 
             * n:          Number of elements in the arrays.
             * I, V, U, Q: Arrays to initialize to zero.
             *
             * withPolarization: If false, only initializes I.
             *                   Otherwise, all arrays are
             *                   initialized to zero.
             */
            template<bool withPolarization>
            void ResetSpectra(
                unsigned int n,
                slibreal_t *I, slibreal_t *Q, slibreal_t *U, slibreal_t *V
            ) {
                unsigned int i;
                for (i = 0; i < n; i++) {
                    I[i] = 0.0;
                    if (withPolarization) {
                        Q[i] = U[i] = V[i] = 0.0;
                    }
                }
            }

            // Pure virtual functions
            virtual slibreal_t Integrate(
                RadiationParticle*, bool,
                slibreal_t*, slibreal_t*,
                slibreal_t*, slibreal_t*
            ) = 0;
    };
}

#endif/*_ANGULAR_DISTRIBUTION_QUADRATURE_2D_H*/
