#ifndef _TEST_ANGULAR_DISTRIBUTION_H
#define _TEST_ANGULAR_DISTRIBUTION_H

#include <softlib/config.h>
#include "UnitTest.h"

class Test_AngularDistribution : public UnitTest {
    private:
    public:
        Test_AngularDistribution(const string& msg) : UnitTest(msg) { }

        bool Run();
};

#endif/*_TEST_ANGULAR_DISTRIBUTION_H*/
