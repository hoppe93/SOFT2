#ifndef _TEST_CONEPROJECTION_REVERSE_H
#define _TEST_CONEPROJECTION_REVERSE_H

#include <string>
#include <softlib/config.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/Cone/Projection/Reverse.h"
#include "../../../../../UnitTest.h"

class Test_ConeProjectionReverse : public UnitTest {
    private:
	public:
		Test_ConeProjectionReverse(const string& msg) : UnitTest(msg) {}

		bool Run(bool);
};

#endif/*_TEST_CONEPROJECTION_REVERSE_H*/
