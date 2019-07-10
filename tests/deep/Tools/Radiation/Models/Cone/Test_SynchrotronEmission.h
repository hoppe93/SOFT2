#ifndef _TEST_SYNCHROTRON_EMISSION_H
#define _TEST_SYNCHROTRON_EMISSION_H

#include <softlib/config.h>
#include "Tools/Radiation/RadiationParticle.h"
#include "../../../../UnitTest.h"

class Test_SynchrotronEmission : public UnitTest {
    private:
        static const unsigned int NTESTPARTICLES;
        static const slibreal_t TESTPARTICLES[][4];
        
        const slibreal_t ANGDIST_TOL = 10.0*SQRT_REAL_EPSILON;
    public:
        Test_SynchrotronEmission(const string& msg) : UnitTest(msg) { }

        bool CheckTotalEmission(const slibreal_t);
        bool CheckSpectrumEmission(const slibreal_t);
        slibreal_t Larmor(__Radiation::RadiationParticle*);
        __Radiation::Detector *GetDetector(unsigned int, slibreal_t l0=4e-7, slibreal_t l1=1e-6);
        slibreal_t GetLambdaC(__Radiation::RadiationParticle*);
        __Radiation::RadiationParticle *GetRadiationParticle(unsigned int, __Radiation::Detector*, MagneticFieldAnalytical2D*);
        bool Run(bool);
};

#endif/*_TEST_SYNCHROTRON_EMISSION_H*/
