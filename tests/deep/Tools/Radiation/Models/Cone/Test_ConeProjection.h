#ifndef _TEST_CONEPROJECTION_H
#define _TEST_CONEPROJECTION_H

#include <string>
#include <softlib/config.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/Cone/Projection/ConeProjection.h"
#include "../../../../UnitTest.h"

struct conemodel_spec {
    slibreal_t x[3];
    slibreal_t vhat[3];
    slibreal_t ppar;
    slibreal_t pperp;
    slibreal_t fraction;
};

class Test_ConeProjection : public UnitTest {
    private:
        static const unsigned int nconespec;
        static const unsigned int conespec_N[];
        static const struct conemodel_spec *conespec[];
        static const slibreal_t *detdir[];
        static const slibreal_t *detpos[];
        static const slibreal_t visang[];
        static const slibreal_t aperture[];
	public:
		Test_ConeProjection(const string& msg) : UnitTest(msg) {}

        bool CompareModelToTable(__Radiation::ConeProjection*, unsigned int, slibreal_t, string&);
        __Radiation::RadiationParticle *GetRadiationParticle(const struct conemodel_spec*, unsigned int, __Radiation::Detector*);
        template<class T>
        bool CompareToAllTables(std::string, slibreal_t);
		bool Run(bool);
        bool RunSOFTv1Comparison();
};

#endif/*_TEST_CONEPROJECTION_H*/
