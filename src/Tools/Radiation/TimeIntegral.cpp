/**
 * Evaluate the time integral.
 */

#include <softlib/config.h>
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

using namespace __Radiation;


void Radiation::HandleTimeIntegral(
    Orbit *o, Particle *p,
    void (Radiation::*evalToroidalIntegral)(
        RadiationParticle&, orbit_type_t,
        slibreal_t, slibreal_t, slibreal_t,
        slibreal_t, slibreal_t
    )
) {
    unsigned int ntau = o->GetNTau();
    slibreal_t *X = o->GetX();
    slibreal_t *P = o->GetP();
    
    unsigned int i, tindex;
    for (i = tindex = 0; i < ntau-1; i++, tindex += 3) {
        RadiationParticle rp(o, i, detector->GetPosition());
        rp.SetDifferentialElements(
            this->dphi, p->GetDRho(),
            p->GetJMomentum1(), p->GetJMomentum2(),
            p->GetJZeta()
        );
        rp.SetDistributionValue(p->GetF());

        model->InitializeTimestep(&rp);

        slibreal_t
            x0 = X[tindex+0],
            y0 = X[tindex+1],
            z  = X[tindex+2],
            px0 = P[tindex+0],
            py0 = P[tindex+1];

        (this->*evalToroidalIntegral)(rp, o->GetOrbitType(), x0, y0, z, px0, py0);
    }
}

