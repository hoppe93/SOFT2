/**
 * Implementation of the unit test for the
 * 'Particle' class in SOFT.
 *
 */

#include <cmath>
#include <cstdlib>
#include <ctime>

// From SOFT
#include "PhaseSpace/Particle.h"
#include "SOFTException.h"

// From SOFT tests
#include "Test_Particle.h"

using namespace std;

const int TEST_PARTICLE_NCOORDS=6;
const int test_particle_coordinates[TEST_PARTICLE_NCOORDS] = {
	Particle::COORDINATE_GAMMA,
	Particle::COORDINATE_P,
	Particle::COORDINATE_PPAR,
	Particle::COORDINATE_PPERP,
	Particle::COORDINATE_THETAP,
	Particle::COORDINATE_XI
};

/**
 * Main routine for testing Particle object.
 * Checks if 'InitializeMomentum()' and its child
 * methods work as expected.
 * 
 * c1:           First set of coordinates; actually tested.
 * c2:			 Second set of coordinates; used for testing area_element.
 * area_element: Area element, p^2 * sin(thetap) * dp * dthetap.
 *
 * RETURNS true if all checks passed. false
 * if one test failed. Prints an error
 * message and returns immediately on
 * failure.
 */
bool Test_Particle::CheckParticleGeneration(
	struct test_particle_coords& c1, struct test_particle_coords& c2,
	struct test_particle_coords& c3
) {
	int i, j, t1, t2;
	slibreal_t p1, p2, dx1dp, dx1dt, dx2dp, dx2dt, dp, dthetap, soft_jacobian, jacobian, ijacobian;
	Particle p;

	for (i = 0; i < TEST_PARTICLE_NCOORDS; i++) {
		for (j = 0; j < TEST_PARTICLE_NCOORDS; j++) {
			if (i == j) continue;

			ResetParticle(p);

			t1 = test_particle_coordinates[i];
			t2 = test_particle_coordinates[j];

			p1 = GetValueFromCoords(t1, c1);
			p2 = GetValueFromCoords(t2, c1);
			
            /* To calculate the Jacobian, we first calculate
             *
             *   dx1/dp  &  dx2/dp,
             *
             * and
             *
             *   dx1/dthetap  &  dx2/dthetap.
             *
             * Using this, we can then compute the Jacobian
             * determinant J that transforms from (p,thetap)
             * to (x1,x2). The inverse of J will then instead
             * transform from (x1,x2) to (p,thetap). And since
             *
             *   d^2p = p^2\sin\thetap dp d\thetap,
             *
             * and
             *
             *   dp d\thetap = (1/J) dx1 dx2,
             *
             * we find
             * 
             *   d^2p = p^2\sin\theta / J dx1 dx2.
             */
			dp = fabs(c2.p - c1.p);
			dthetap = fabs(c3.thetap - c1.thetap);

			dx1dp = (GetValueFromCoords(t1, c2) - p1) / dp;
			dx1dt = (GetValueFromCoords(t1, c3) - p1) / dthetap;

			dx2dp = (GetValueFromCoords(t2, c2) - p2) / dp;
			dx2dt = (GetValueFromCoords(t2, c3) - p2) / dthetap;

			// We calculate 1/J first and then invert the result
			ijacobian = fabs((dx1dp*dx2dt) - (dx2dp*dx1dt));
			jacobian = c1.p*c1.pperp / ijacobian;

			try {
				p.InitializeMomentum(t1, t2, p1, p2, 1.0, 1.0);
				soft_jacobian = p.GetJMomentum1() * p.GetJMomentum2();

				if (fabs((c1.gamma - p.GetGamma())/c1.gamma) > TEST_PARTICLE_LIMIT) {
					this->PrintError(
						"Coordinate gamma: Results don't match. Input: %s / %s. Delta = %e",
						Particle::GetCoordinateName(t1),
						Particle::GetCoordinateName(t2),
						fabs((c1.gamma - p.GetGamma())/c1.gamma)
					);
					return false;
				} else if (fabs((c1.p - p.GetMomentum())/c1.p) > TEST_PARTICLE_LIMIT) {
					this->PrintError(
						"Coordinate p: Results don't match. Input: %s / %s. Delta = %e",
						Particle::GetCoordinateName(t1),
						Particle::GetCoordinateName(t2),
						fabs((c1.p - p.GetMomentum())/c1.p)
					);
					return false;
				} else if (fabs((c1.ppar - p.GetPpar())/c1.ppar) > TEST_PARTICLE_LIMIT) {
					this->PrintError(
						"Coordinate ppar: Results don't match. Input: %s / %s. Delta = %e",
						Particle::GetCoordinateName(t1),
						Particle::GetCoordinateName(t2),
						fabs((c1.ppar - p.GetPpar())/c1.ppar)
					);
					return false;
				} else if (fabs((c1.pperp - p.GetPperp())/c1.pperp) > TEST_PARTICLE_LIMIT) {
					this->PrintError(
						"Coordinate pperp: Results don't match. Input: %s / %s. Delta = %e. par2 = %e",
						Particle::GetCoordinateName(t1),
						Particle::GetCoordinateName(t2),
						fabs((c1.pperp - p.GetPperp())/c1.pperp),
						p2
					);
					return false;
				} else if (fabs((c1.thetap - p.GetThetap())/c1.thetap) > TEST_PARTICLE_LIMIT) {
					this->PrintError(
						"Coordinate thetap: Results don't match. Input: %s / %s. Delta = %e",
						Particle::GetCoordinateName(t1),
						Particle::GetCoordinateName(t2),
						fabs((c1.thetap - p.GetThetap())/c1.thetap)
					);
					return false;
				} else if (fabs((c1.xi - p.GetXi())/c1.xi) > TEST_PARTICLE_LIMIT) {
					this->PrintError(
						"Coordinate xi: Results don't match. Input: %s / %s. Delta = %e",
						Particle::GetCoordinateName(t1),
						Particle::GetCoordinateName(t2),
						fabs((c1.xi - p.GetXi())/c1.xi)
					);
					return false;
				} else if (fabs((jacobian - soft_jacobian)/jacobian) > TEST_PARTICLE_LIMIT_JAC) {
					this->PrintError(
						"Jacobian: Results don't match. Input: %s / %s. Delta = %e.",
						Particle::GetCoordinateName(t1),
						Particle::GetCoordinateName(t2),
						fabs((jacobian - soft_jacobian)/jacobian)
					);
					return false;
				}
			} catch (SOFTException& e) {
				// This means the coordinate combination was
				// unsupported, which is fine and can be ignored.
			}
		}
	}

	return true;
}

/**
 * Calculate xi, ppar, pperp and gamma
 * from a 'coords' object whose p and thetap
 * components have been initialized.
 *
 * c: test_particle_coords object to populate.
 */
void Test_Particle::GenerateCoordinates(struct test_particle_coords& c) {
	c.xi    = cos(c.thetap);
	c.ppar  = c.p * c.xi;
	c.pperp = c.p * sin(c.thetap);
	c.gamma = hypot(c.p, 1.0);
}

/**
 * Returns the entry of a 'struct test_particle_coords'
 * that corresponds to the given Particle::COORDINATE_???
 * coordinate type.
 *
 * t: Coordinate type to return.
 * c: Coordinate object from which to pick a value.
 */
slibreal_t Test_Particle::GetValueFromCoords(const int t, const struct test_particle_coords& c) {
	switch (t) {
		case Particle::COORDINATE_GAMMA: return c.gamma;
		case Particle::COORDINATE_P: return c.p;
		case Particle::COORDINATE_PPAR: return c.ppar;
		case Particle::COORDINATE_PPERP: return c.pperp;
		case Particle::COORDINATE_THETAP: return c.thetap;
		case Particle::COORDINATE_XI: return c.xi;
		default: throw SOFTLibException("Unrecognized coordinate type: "+to_string(t));
	}
}

/**
 * Reset a SOFT particle, setting all its
 * momentum properties to NaN.
 *
 * p: Particle object to reset.
 */
void Test_Particle::ResetParticle(Particle& p) {
	p.SetGamma(nan(""));
	p.SetMomentum(nan(""));
	p.SetPpar(nan(""));
	p.SetPperp(nan(""));
	p.SetThetap(nan(""));
	p.SetXi(nan(""));
	p.SetDMomentum1(nan(""));
	p.SetDMomentum2(nan(""));
}
/**
 * Run the unit test.
 */
bool Test_Particle::Run(bool) {
	int i, failed;
	struct test_particle_coords c1, c2, c3;
	slibreal_t dp, dthetap;

	dp = TEST_PARTICLE_DP;
	dthetap = TEST_PARTICLE_DTHETAP;

	failed = 0;
	for (i = 0; i < TEST_PARTICLE_NTESTS; i++) {
		// Generate test point
		c1.p = TEST_PARTICLE_PMAX * Rand();
		c1.thetap = M_PI * Rand();

		c2.p = c1.p + dp;
		c2.thetap = c1.thetap;
		c3.p = c1.p;
		c3.thetap = c1.thetap + dthetap;

		// If thetap > pi, redo step
		if (c3.thetap > M_PI) {
			i--;
			continue;
		}

		GenerateCoordinates(c1);
		GenerateCoordinates(c2);
		GenerateCoordinates(c3);

		if (!CheckParticleGeneration(c1, c2, c3))
			failed++;
	}

	return (failed < TEST_PARTICLE_NTESTS*0.01);
}

