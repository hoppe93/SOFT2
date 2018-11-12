#ifndef _TEST_ANGULAR_DISTRIBUTION_DRIFTS_H
#define _TEST_ANGULAR_DISTRIBUTION_DRIFTS_H

#include <softlib/config.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include "Tools/Radiation/RadiationParticle.h"
#include "GOAMomentum.h"
#include "../../../../UnitTest.h"

class Test_AngularDistributionDrifts : public UnitTest {
    struct betadot_params {
        GOAMomentum *goa;
        slibreal_t ppar, pperp;
    };

    private:
        static const unsigned int NTESTPARTICLES;
        static const slibreal_t TESTPARTICLES[][3];
        
        const slibreal_t ANGDIST_TOL = 10.0*SQRT_REAL_EPSILON;

        static double betaDotSquared(double, void*);
        struct intpar {
            slibreal_t betapar;
            __Radiation::ADSynchrotronEmission *ade;
        };
    public:
        Test_AngularDistributionDrifts(const string& msg) : UnitTest(msg) { }

        bool Run(bool);

        bool CheckAngularDistributionIntegral(
            __Radiation::ADSynchrotronEmission&, __Radiation::RadiationParticle*,
            MagneticFieldAnalytical2D*, slibreal_t
        );
        slibreal_t EvaluateMuIntegral(__Radiation::ADSynchrotronEmission&, __Radiation::RadiationParticle*, bool adapt=true);
        static double EvaluateMuIntegral_inner(double, void*);
        __Radiation::Detector *GetDetector(unsigned int, slibreal_t l0=4e-7, slibreal_t l1=1e-6);
        Vector<3> GetGuidingCenterMomentum(
            MagneticField2D*, struct magnetic_field_data&,
            const Vector<3>&, const slibreal_t, const slibreal_t
        );
        __Radiation::RadiationParticle *GetRadiationParticle(MagneticFieldAnalytical2D*, __Radiation::Detector*, unsigned int);
        slibreal_t Larmor(__Radiation::RadiationParticle*, MagneticFieldAnalytical2D*);

        bool CompareBetas(const slibreal_t);
        bool CompareBetaDots(const slibreal_t);
        bool VerifyIntegral(const slibreal_t);
        bool VerifyP2Scaling(const slibreal_t);
};

#endif/*_TEST_ANGULAR_DISTRIBUTION_DRIFTS_H*/
