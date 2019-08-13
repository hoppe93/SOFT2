#ifndef _TEST_BREMSSTRHLUNG_SCREENED_EMISSION_H
#define _TEST_BREMSSTRAHLUNG_SCREENED_EMISSION_H

#include <softlib/config.h>
#include "Tools/Radiation/RadiationParticle.h"
#include "../../../../../UnitTest.h"

class Test_BremsstrahlungScreenedEmission : public UnitTest {
    private:
        static const unsigned int NTESTVALUES;
        static const slibreal_t TEST_Z[][6];
        static const slibreal_t TEST_Z0[][6];
        static const slibreal_t TEST_DENS[][6];
        static const slibreal_t TEST_GAMMA_INTSPEC[][2];
        const slibreal_t ANGDIST_TOL = 10.0*SQRT_REAL_EPSILON;
    public:
        Test_BremsstrahlungScreenedEmission(const string& msg) : UnitTest(msg) { }

        bool CheckTotalEmission(const slibreal_t);
        bool CheckSpectrumEmission(const slibreal_t);
        __Radiation::Detector *GetDetector(unsigned int, slibreal_t l0=4e-7, slibreal_t l1=1e-6);
        __Radiation::RadiationParticle *GetRadiationParticle(unsigned int, __Radiation::Detector*, MagneticFieldAnalytical2D*);
        bool Run(bool);
};

#endif/*_TEST_BREMSSTRAHLUNG_SCREENED_EMISSION_H*/
