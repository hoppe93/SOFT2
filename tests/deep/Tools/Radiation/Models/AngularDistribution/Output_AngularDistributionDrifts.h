#ifndef _OUTPUT_ANGULAR_DISTRIBUTION_DRIFTS_H
#define _OUTPUT_ANGULAR_DISTRIBUTION_DRIFTS_H

#include <softlib/config.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "../../../../UnitTest.h"

class Output_AngularDistributionDrifts : public UnitTest {
    private:
        static const unsigned int NTESTPARTICLES;
        static const slibreal_t TESTPARTICLES[][3];
    public:
        Output_AngularDistributionDrifts(const string& msg) : UnitTest(msg) { }

        __Radiation::Detector *GetDetector(unsigned int, slibreal_t l0=400e-9, slibreal_t l1=1000e-9);
        __Radiation::RadiationParticle *GetRadiationParticle(MagneticFieldAnalytical2D*, __Radiation::Detector*, unsigned int);
        Vector<3> GetGuidingCenterMomentum(MagneticField2D*, struct magnetic_field_data&, const Vector<3>&, const slibreal_t, const slibreal_t);
        void OutputIntegrand();
        bool Run(bool);
};

#endif/*_OUTPUT_ANGULAR_DISTRIBUTION_DRIFTS_H*/
