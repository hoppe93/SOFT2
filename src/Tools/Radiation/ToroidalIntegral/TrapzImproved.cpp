/**
 * Handle the toroidal integral using a regular
 * trapezoidal rule.
 */

#include <cstring>
#include <softlib/config.h>
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

using namespace __Radiation;
using namespace std;

void Radiation::HandleTrapzImproved(Orbit *o, Particle *p) {
    unsigned int ntau, i, tindex;
    slibreal_t z, x0, y0, px0, py0, *X, *P;

    ntau = o->GetNTau();
    X = o->GetX();
    P = o->GetP();

    model->InitializeOrbit(o);

    // tau integral
    // (ntau-1, since the first and last points are the same)
    for (i = tindex = 0; i < ntau-1; i++, tindex += 3) {
        RadiationParticle rp(o, i, detector->GetPosition());
        rp.SetDifferentialElements(
            this->dphi, p->GetDRho(),
            p->GetJMomentum1(), p->GetJMomentum2()
        );
        rp.SetDistributionValue(p->GetF());

        model->InitializeTimestep(&rp);

        x0  = X[tindex+0];
        y0  = X[tindex+1];
        z   = X[tindex+2];

        px0 = P[tindex+0];
        py0 = P[tindex+1];

        memset(this->torflags, 0, ntoroidal);

        unsigned int phi1, phi2;
        LocateSurfaceOfVisibility(&rp, &phi1, &phi2);

        slibreal_t mx = 0;
        IntegrateToroidalImproved(rp, x0, y0, z, px0, py0, phi1, -1, mx);
        IntegrateToroidalImproved(rp, x0, y0, z, px0, py0, phi1+1, +1, mx);

        if (phi2 != phi1) {
            IntegrateToroidalImproved(rp, x0, y0, z, px0, py0, phi2, -1, mx);
            IntegrateToroidalImproved(rp, x0, y0, z, px0, py0, phi2+1, +1, mx);
        }
    }
}

unsigned int Radiation::IntegrateToroidalImproved(
    RadiationParticle &rp, slibreal_t x0, slibreal_t y0, slibreal_t z,
    slibreal_t px0, slibreal_t py0, int startj, int sgn, slibreal_t &mx
) {
    slibreal_t x, y, px, py;
    bool has_outer_wall = (wall_opacity!=WALL_OPACITY_SEMI_TRANSPARENT);
    Vector<3> rcp, detx;
    unsigned int niter = 0;

    detx = detector->GetPosition();

    // Toroidal integral
    int j = startj;
    if (j == (int)ntoroidal)
        j = 0;

    do {
        // Has this point already been counted?
        if (this->torflags[j])
            goto NEXT;
        else this->torflags[j] = 1;

        x = x0*cosphi[j] + y0*sinphi[j];
        y =-x0*sinphi[j] + y0*cosphi[j];

        px = px0*cosphi[j] + py0*sinphi[j];
        py =-px0*sinphi[j] + py0*cosphi[j];

        // Check if within FOV
        if (!IsWithinFieldOfView(x, y, z, rcp)) {
            if (j != startj)
                break;
            goto NEXT;
        }

        rp.UpdateXY(x, y, px, py, rcp);

        model->HandleParticle(&rp, sinphi[j], cosphi[j]);

        if (!model->IsNonZero() || model->GetPower() < mx*this->torthreshold) {
            break;
        } else if (wall_opacity == WALL_OPACITY_TRANSPARENT ||
            !magfield->IntersectsDomain3D(
            x, y, z, detx[0], detx[1], detx[2], has_outer_wall)
        ) {
            if (model->GetPower() > mx)
                mx = model->GetPower();
            RegisterOutput(&rp);
        }

NEXT:
        j += sgn;
        if (sgn > 0 && j >= (int)ntoroidal)
            j = 0;
        else if (sgn < 0 && j < 0)
            j = ntoroidal-1;

        niter++;
    } while (j != startj);

    return niter;
}

/**
 * Reset toroidal flags array.
 */
void Radiation::ResetToroidalFlags() {
    memset(this->torflags, 0, ntoroidal*sizeof(char));
}

