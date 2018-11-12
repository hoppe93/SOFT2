#ifndef _TEST_ISOTROPIC_H
#define _TEST_ISOTROPIC_H

#include <softlib/config.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/Isotropic.h"
#include "../../../../UnitTest.h"

class Test_Isotropic : public UnitTest {
    private:
        // Simpson's rule is relatively slow to converge,
        // so we will need a rather high tolerance
        const slibreal_t ISO_TOL = 5e-3;
        const unsigned int NSIMPSON_STEPS = 10;
	public:
		Test_Isotropic(const string& msg) : UnitTest(msg) {}

        bool CheckIsotropicEmission(__Radiation::Isotropic*);
        slibreal_t NumericIntegral(__Radiation::Detector*, __Radiation::RadiationParticle*);
        slibreal_t NumericIntegral_inner(
            slibreal_t, slibreal_t, slibreal_t,
            Vector<3>&, Vector<3>&, Vector<3>&, Vector<3>&
        );
        slibreal_t NumericIntegral_intg(
            slibreal_t, slibreal_t,
            Vector<3>&, Vector<3>&, Vector<3>&, Vector<3>&
        );
		bool Run(bool);
};

#endif/*_TEST_ISOTROPIC_H*/
