#ifndef _TEST_PARTICLE_PUSHER_H
#define _TEST_PARTICLE_PUSHER_H

#include <softlib/config.h>
#include "Orbit/Orbit.h"
#include "../../UnitTest.h"

#define NORBITPOINTS 100
struct orbit {
    slibreal_t r, p, thetap;
    unsigned int elems;
    slibreal_t *t, *x, *y, *z;
};

class Test_ParticlePusher : public UnitTest {
    private:
        slibreal_t
            ORBIT_TOLERANCE = 1e-12,
            GCMA_TOLERANCE=100.0*ORBIT_TOLERANCE,
            GC_TOLP=1000.0*ORBIT_TOLERANCE,
            GC_TOLMU=10.0*ORBIT_TOLERANCE,
            GC_SOL_TOL=1e-8,
            PART_SOL_TOL=1e-8;


        static const unsigned int ngcorbits, ngcorbits_nd, npartorbits,
            guiding_center_orbits_nmax, guiding_center_orbits_nd_nmax, particle_orbits_nmax;
        static struct orbit guiding_center_orbits[];
        static struct orbit guiding_center_orbits_nd[];
        static struct orbit particle_orbits[];

        const slibreal_t
            B0=5.0,
            Rm=0.68,
            rminor=0.22,
            zaxis=0;

    public:
        Test_ParticlePusher(const string& name) : UnitTest(name) {}

        bool Run(bool);

        bool GCMagneticAxis();
        bool GuidingCenterEquation();
        bool GuidingCenterEquation_Conservation(bool);
        bool GuidingCenterEquation_TableCompare(bool);
        Orbit *GenerateOrbit(const slibreal_t, const slibreal_t, const slibreal_t, bool drifts=false, const string& eqn="guiding-center");
        ParticlePusher *GenerateOrbit_push(const slibreal_t, const slibreal_t, const slibreal_t, bool drifts=false, const string& eqn="guiding-center");
};

#endif/*_TEST_PARTICLE_PUSHER_H*/
