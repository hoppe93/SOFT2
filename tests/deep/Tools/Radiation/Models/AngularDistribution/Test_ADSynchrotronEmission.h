#ifndef _TEST_ADSYNCHROTRON_EMISSION_H
#define _TEST_ADSYNCHROTRON_EMISSION_H

#include <softlib/config.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"
#include "../../../../UnitTest.h"

class Test_ADSynchrotronEmission : public UnitTest {
    private:
        static const unsigned int NTESTPARTICLES;
        static const slibreal_t TESTPARTICLES[][3];
        
        const slibreal_t ANGDIST_TOL = 10.0*SQRT_REAL_EPSILON;
    public:
        Test_ADSynchrotronEmission(const std::string& msg) : UnitTest(msg) { }

        bool CheckAngularDistribution(
            __Radiation::ADSynchrotronEmission&,
            __Radiation::RadiationParticle*,
            slibreal_t
        );
        bool CheckAngularSpectralDistribution(
            __Radiation::ADSynchrotronEmission&,
            __Radiation::RadiationParticle*,
            slibreal_t
        );
        slibreal_t GetLambdaC(__Radiation::RadiationParticle*);
        slibreal_t Larmor(__Radiation::RadiationParticle*);
        __Radiation::Detector *GetDetector(unsigned int, slibreal_t l0=4e-7, slibreal_t l1=1e-6);
        __Radiation::RadiationParticle *GetRadiationParticle(unsigned int, __Radiation::Detector*);
        bool Run(bool);
};

#endif/*_TEST_ADSYNCHROTRON_EMISSION_H*/
