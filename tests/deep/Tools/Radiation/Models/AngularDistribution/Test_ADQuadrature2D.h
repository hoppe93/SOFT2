#ifndef _TEST_ANGULAR_DISTRIBUTION_QUADRATURE2D_H
#define _TEST_ANGULAR_DISTRIBUTION_QUADRATURE2D_H

#include <softlib/config.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/AngularDistribution/ADQuadrature2D.h"
#include "../../../../UnitTest.h"

class Test_ADQuadrature2D_Emission : public __Radiation::ADEmission {
    private:
    public:
        Test_ADQuadrature2D_Emission(__Radiation::Detector *d) : ADEmission(d) {}

        virtual slibreal_t Evaluate(Vector<3>&, slibreal_t, slibreal_t, bool) override;
        virtual void InitializeToroidalStep(const slibreal_t, const slibreal_t) override {}
        virtual void Prepare(__Radiation::RadiationParticle*, bool) override {}

        virtual void CalculateAngularDistribution(Vector<3>&, slibreal_t, slibreal_t) override;
        virtual void CalculatePolarization(Vector<3>&, slibreal_t, slibreal_t) override;
        virtual void CalculateSpectrum(Vector<3>&, slibreal_t, slibreal_t) override;
};

class Test_ADQuadrature2D : public UnitTest {
    private:
        static const slibreal_t DETECTOR_APERTURE;
        static const unsigned int NTESTPARTICLES;
        static const slibreal_t test_particles[][3];
    public:
        Test_ADQuadrature2D(const std::string& msg) : UnitTest(msg) { }

        __Radiation::Detector *GetDetector(unsigned int);
        __Radiation::RadiationParticle *GetRadiationParticle(unsigned int, __Radiation::Detector*);
        template<class T, typename ... Args> bool IsotropicEmission(slibreal_t, Args&& ...);
        bool IsotropicEmissionTest(__Radiation::ADQuadrature2D&, __Radiation::RadiationParticle*, slibreal_t);
        bool IsotropicEmissionVerify(
            __Radiation::ADQuadrature2D&, __Radiation::RadiationParticle*, bool, slibreal_t,
            slibreal_t*, slibreal_t*, slibreal_t*, slibreal_t*
        );
        bool Run(bool);
};

#endif/*_TEST_ANGULAR_DISTRIBUTION_QUADRATURE2D_H*/
