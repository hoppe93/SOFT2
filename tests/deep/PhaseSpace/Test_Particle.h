#ifndef _TEST_PARTICLE_H
#define _TEST_PARTICLE_H

#include <softlib/config.h>
#include "PhaseSpace/Particle.h"
#include "../UnitTest.h"

#define TEST_PARTICLE_DP        sqrt(REAL_EPSILON)	// Distance in p between the two test points
#define TEST_PARTICLE_DTHETAP   sqrt(REAL_EPSILON)	// Distance in thetap between the two test points
#define TEST_PARTICLE_LIMIT     sqrt(REAL_EPSILON)  // Limit below which two numbers should be considered equal
#define TEST_PARTICLE_LIMIT_JAC 1e-4                // Limit on the jacobian
#define TEST_PARTICLE_NTESTS    1000				// Number of test iterations
#define TEST_PARTICLE_PMAX      100.0				// Maximum allowed value of p (momentum, normalized to mc)

struct test_particle_coords { slibreal_t gamma, p, ppar, pperp, thetap, xi; };

class Test_Particle : public UnitTest {
	public:
		Test_Particle(const string& msg) : UnitTest(msg) {}

		bool CheckParticleGeneration(struct test_particle_coords&, struct test_particle_coords&, struct test_particle_coords&);
		void GenerateCoordinates(struct test_particle_coords&);
		slibreal_t GetValueFromCoords(const int, const struct test_particle_coords&);
		void ResetParticle(Particle&);
		bool Run(bool);
};

#endif/*_TEST_PARTICLE_H*/
