#ifndef _ANGULAR_DISTRIBUTION_SYNCHROTRON_EMISSION_H
#define _ANGULAR_DISTRIBUTION_SYNCHROTRON_EMISSION_H

#include <gsl/gsl_integration.h>
#include <softlib/config.h>
#include <softlib/Vector.h>
#include "SOFT.h"
#include "Tools/Radiation/Models/AngularDistributionException.h"
#include "Tools/Radiation/Models/AngularDistribution/ADEmission.h"
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationException.h"
#include "Tools/Radiation/RadiationParticle.h"

slibreal_t synchrotron_func(const slibreal_t);

namespace __Radiation {
    class ADSynchrotronEmission : public ADEmission {
        public:
        struct angdist_params {
            Vector<3> n, Bvec, bHat, oneHat, twoHat,
                PGC, curlBhat, curlB, gradB,
                oneDotNablaB, twoDotNablaB;

            // Particle velocity coefficients
            Vector<3> sin, cos, sincos, sin2, cos2;
            slibreal_t B, ppar, pperp, gamma, m, q,
                iGamma, betaDotFactor, rho, gammapar2;
            slibreal_t sinphi, cosphi;      // Toroidal angle
        };

        private:
            bool includeDrifts = false;

            static const slibreal_t ANGDIST_PREFAC_ND;
            
            double qagsEpsAbs=0.0,      // Absolute tolerance currently disabled
                   qagsEpsRel=1e-3;
            std::size_t qagsLimit = 100;
            gsl_integration_workspace *qagsWorkspace=nullptr;

            // Internal cache variables
            slibreal_t 
                B, beta, gamma, gamma3, igamma, igamma2, sinThetap, cosThetap,
                p, pperp, prefactor, lambdac, betapar;
            struct angdist_params angdistParams;

            static void _Rotate(Vector<3>&, const slibreal_t, const slibreal_t);
        public:
            ADSynchrotronEmission(
                Detector*, struct global_settings*,
                std::size_t qagsLimit=100, slibreal_t qagsEpsRel=1e-3
            );
            ~ADSynchrotronEmission();

            virtual slibreal_t Evaluate(Vector<3>&, slibreal_t, slibreal_t, bool) override;

            virtual void CalculateAngularDistribution(Vector<3>&, slibreal_t, slibreal_t) override;
            virtual void CalculatePolarization(Vector<3>&, slibreal_t, slibreal_t) override;
            virtual void CalculateSpectrum(Vector<3>&, slibreal_t, slibreal_t) override;
            template<bool calculatePolarization>
            void __CalculateSpectrum(Vector<3>&, slibreal_t, slibreal_t);

            struct angdist_params *GetParams() { return &angdistParams; }

            virtual void InitializeToroidalStep(const slibreal_t, const slibreal_t) override;
            void InitializeToroidalStepAD(const slibreal_t, const slibreal_t);
            void IntegrateSpectrum();

            virtual void Prepare(RadiationParticle*, bool) override;
            void PrepareAngularDistribution(RadiationParticle*);
            void PreparePolarization(RadiationParticle*);
            void PrepareSpectrum(RadiationParticle*);

            void RotateOneTwo(const slibreal_t);

            slibreal_t CalculateAngularDistribution_ZerothOrderGC(slibreal_t, slibreal_t);
            slibreal_t CalculateAngularDistribution_FirstOrderGC(Vector<3>&, slibreal_t, slibreal_t);
            static void PrepareFirstOrder(RadiationParticle*, struct angdist_params*, slibreal_t*);

            // Internal routines for 1st order GC theory
            static Vector<3> &CalculateBeta(struct angdist_params*, slibreal_t, slibreal_t, Vector<3>&);
            static Vector<3> &CalculatePGC(struct angdist_params*);
            static Vector<3> &CalculateBetaDot(struct angdist_params*, slibreal_t, slibreal_t, const Vector<3>&, Vector<3>&);

            static double synchrotron_first_order_gc(double, void*);
    };

    class ADSynchrotronException : public AngularDistributionException {
        public:
            template<typename ... Args>
            ADSynchrotronException(const std::string &msg, Args&& ... args)
                : AngularDistributionException(msg, std::forward<Args>(args) ...) {
                AddModule("Synchrotron");
            }
    };
}

#endif/*_ANGULAR_DISTRIBUTION_SYNCHROTRON_EMISSION_H*/
