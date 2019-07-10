/* Unit tests */

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <softlib/config.h>

#include "UnitTest.h"

// Tests
#include "Orbit/ParticlePusher/Test_ParticlePusher.h"
#include "PhaseSpace/Test_Particle.h"
#include "PhaseSpace/Test_ParticleGenerator.h"
#include "Tools/Radiation/Models/AngularDistribution/Output_AngularDistributionDrifts.h"
#include "Tools/Radiation/Models/AngularDistribution/Test_ADQuadrature2D.h"
#include "Tools/Radiation/Models/AngularDistribution/Test_ADSynchrotronEmission.h"
#include "Tools/Radiation/Models/AngularDistribution/Test_AngularDistributionDrifts.h"
#include "Tools/Radiation/Models/Cone/Test_ConeProjection.h"
#include "Tools/Radiation/Models/Cone/Test_SynchrotronEmission.h"
#include "Tools/Radiation/Models/Cone/Bremsstrahlung/Test_BremsstrahlungScreenedEmission.h" // new bremsstrahlung test
#include "Tools/Radiation/Models/Cone/Reverse/Test_ConeProjectionReverse.h"
#include "Tools/Radiation/Models/Isotropic/Test_Isotropic.h"


using namespace std;

vector<UnitTest*> tests;

void add_test(UnitTest *t) {
	tests.push_back(t);
}
void init() {
	add_test(new Test_Particle("particle"));
    add_test(new Test_ParticleGenerator("particle_generator"));
    add_test(new Test_ParticlePusher("particle_pusher"));
    add_test(new Test_ADQuadrature2D("radiation_model_angdist_quad2d"));
    add_test(new Test_ADSynchrotronEmission("radiation_model_angdist_synchrotron"));
    add_test(new Output_AngularDistributionDrifts("radiation_model_angdist_synchrotron_drifts_output"));
    add_test(new Test_AngularDistributionDrifts("radiation_model_angdist_synchrotron_drifts"));
    add_test(new Test_ConeProjection("radiation_model_cone_projection"));
    add_test(new Test_ConeProjectionReverse("radiation_model_cone_projection_reverse"));
    add_test(new Test_SynchrotronEmission("radiation_model_cone_synchrotron"));
    add_test(new Test_BremsstrahlungScreenedEmission("radiation_model_cone_bremsstrahlung_screened")); //new bremsstrahlung test
    add_test(new Test_Isotropic("radiation_model_isotropic"));
}

/**
 * Check if a test with the given name
 * exists in the test list "tests".
 *
 * name: Name of test to look for.
 *
 * RETURNS the index of test if found,
 * otherwise, -1.
 */
int has_test(const char *name) {
	size_t i;
	for (i = 0; i < tests.size(); i++) {
		if (tests[i]->HasName(name))
			return i;
	}

	return -1;
}

/**
 * Run the test with index 'index' in
 * the test list "tests".
 *
 * index: Index in "tests" of the test
 *    to run.
 * runningAll: TRUE if all (or a group of)
 *    tests are being run.
 *
 * RETURNS 1 if the test failed, 0 otherwise.
 */
int run_test(int index, bool runningAll) {
    cout << "\x1B[1m:: " << tests[index]->GetName() << "\x1B[0m" << endl;
    try {
        if (tests[index]->Run(runningAll)) {
            cout << "\x1B[1;32m[SUCCESS]\x1B[0m Test '" << tests[index]->GetName() << "' completed successfully." << endl;
            return 0;
        } else {
            cout << "\x1B[1;31m[FAIL]\x1B[0m    Test '" << tests[index]->GetName() << "' failed." << endl;
            return 1;
        }
    } catch (SOFTLibException& ex) {
        cout << "\x1B[1;31m[ERROR]\x1B[0m   --> " << ex.whats() << endl;
        cout << "\x1B[1;31m[FAIL]\x1B[0m    Test '" << tests[index]->GetName() << "' unexpectedly threw exception." << endl;
    }

    return 1;
}
/**
 * Runs all tests in the list "tests".
 *
 * RETURNS the number of tests that failed.
 */
int run_all_tests() {
	size_t i;
	int failed = 0;

	for (i = 0; i < tests.size(); i++) {
		failed += run_test(i, true);
	}

	return failed;
}

/**
 * Print help text.
 */
void help() {
    printf(
        "Unit tests for SOFT\n\n"

        "Usage:\n"
        "    soft_tests          Show this help message.\n"
        "    soft_tests all      Run all tests.\n"
        "    soft_tests [test1 [test2 [...]]]\n"
        "                        Runs the tests with names 'test1', 'test2' etc.\n\n"

        "Available tests:\n"
    );
    for (unsigned int i = 0; i < tests.size(); i++) {
        printf("    %s\n", tests[i]->GetName().c_str());
    }

    printf("\n");
}

int main(int argc, char *argv[]) {
	int i, t, failed=0, total=0;

	init();

	if (argc == 1) {
        help();
    } else if (argc == 2 && !strcmp(argv[1], "all")) {
		failed = run_all_tests();
        total = tests.size();
	} else {
        total = argc-1;
		for (i = 1, failed=0; i < argc; i++) {
			if ((t=has_test(argv[i]))>=0) {
				failed += run_test(t, (argc>2));
			} else {
				cerr << "\x1B[1;31m[ERROR]\x1B[0m   Unrecognized test: '" << argv[i] << "'." << endl;
				failed++;
			}
		}
	}

	cout << (total-failed) << " of " << total << " tests passed." << endl;

	return failed;
}

