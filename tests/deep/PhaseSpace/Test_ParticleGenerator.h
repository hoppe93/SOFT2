#ifndef _TEST_PARTICLE_GENERATOR_H
#define _TEST_PARTICLE_GENERATOR_H

#include <softlib/config.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include "Orbit/GuidingCenterEquation.h"
#include "PhaseSpace/ParticleGenerator.h"
#include "../UnitTest.h"

#define TEST_PARTICLE_GENERATOR_NRAD        10
#define TEST_PARTICLE_GENERATOR_NP1         10
#define TEST_PARTICLE_GENERATOR_NP2         10
#define TEST_PARTICLE_GENERATOR_PMAX        100

class Test_ParticleGenerator : public UnitTest {
    private:
        static const slibreal_t
            effax_ppar0, effax_ppar1,
            effax_pperp0, effax_pperp1;
        static const unsigned int
            effax_nppar, effax_npperp;
	public:
		Test_ParticleGenerator(const string& msg) : UnitTest(msg) {}

        bool CheckEffectiveMagneticAxis();
        bool CheckPhaseSpaceGrid();
        bool Generate(
            MagneticFieldAnalytical2D*, struct global_settings*,
            const string&, const string&, const string&,
            const slibreal_t, const slibreal_t, const slibreal_t,
            const slibreal_t, const slibreal_t, const slibreal_t
        );
        slibreal_t GetDriftShift(MagneticField2D*, GuidingCenterEquation*, const slibreal_t, const slibreal_t, const slibreal_t);
        slibreal_t GetValueFromParamName(Particle*, const string&, MagneticFieldAnalytical2D*);
        bool TestInvalidInput();
		bool Run(bool);
};

#endif/*_TEST_PARTICLE_GENERATOR_H*/
