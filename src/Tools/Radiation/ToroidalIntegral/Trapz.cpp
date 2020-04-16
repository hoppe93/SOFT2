/**
 * Handle the toroidal integral using a regular
 * trapezoidal rule.
 */

#include <softlib/config.h>
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

using namespace __Radiation;

void Radiation::HandleTrapz(Orbit *o, Particle *p) {
    unsigned int ntau, i, j, tindex;
    slibreal_t x, y, z, x0, y0, z0, px, py, px0, py0, *X, *P;
    Vector<3> rcp, detx;
    bool has_outer_wall = (wall_opacity!=WALL_OPACITY_SEMI_TRANSPARENT);

    ntau = o->GetNTau();
    X = o->GetX();
    P = o->GetP();
    detx = detector->GetPosition();

    model->InitializeOrbit(o);
    orbit_type_t otype = o->GetOrbitType(); //gissning

    // tau integral
    // (ntau-1, since the first and last points are the same)
    for (i = tindex = 0; i < ntau-1; i++, tindex += 3) {
        RadiationParticle rp(o, i, detector->GetPosition());
        rp.SetDifferentialElements(
            this->dphi, p->GetDRho(),
            p->GetJMomentum1(), p->GetJMomentum2(),
            p->GetJZeta()
        );
        rp.SetDistributionValue(p->GetF());

        model->InitializeTimestep(&rp);

        x0  = X[tindex+0];
        y0  = X[tindex+1];
        z0  = X[tindex+2];

        px0 = P[tindex+0];
        py0 = P[tindex+1];

        // Toroidal integral
        for (j = 0; j < ntoroidal; j++) {
            x = x0*cosphi[j] + y0*sinphi[j];
            y =-x0*sinphi[j] + y0*cosphi[j];

            px = px0*cosphi[j] + py0*sinphi[j];
            py =-px0*sinphi[j] + py0*cosphi[j];

            // If drifts are enabled and we're using a guiding-center cone
            // model, then we should also shift the guiding-center position
            // to take the finite Larmor radius into account
            /*if (this->shiftLarmorRadius) {
                z = z0;
                this->ShiftLarmorRadius(x, y, z, px, py, rp.GetP()[2], &rp);
            }*/

            // Check if within FOV
            if (!IsWithinFieldOfView(x, y, z, rcp))
                continue;

            rp.UpdateXY(x, y, px, py, rcp);

            model->HandleParticle(&rp, otype, sinphi[j], cosphi[j]);

            if (!model->IsNonZero()) continue;
            else if (wall_opacity == WALL_OPACITY_TRANSPARENT ||
                !magfield->IntersectsDomain3D(
                x, y, z, detx[0], detx[1], detx[2], has_outer_wall)
            ) {
                RegisterOutput(&rp);
            }
        }
    }
}

/**
 * Evaluate only the toroidal integral.
 *
 * rp:        Radiation particle object
 *            representing the particle state.
 * otype:     Type of orbit passed to this function
 *            (particle or guiding-center).
 * x0, y0, z: Position of the particle at the
 *            initial toroidal angle.
 * px0, py0:  Momentum of the particle at the
 *            initial toroidal angle.
 */
void Radiation::EvaluateToroidalTrapz(
    RadiationParticle &rp, orbit_type_t otype,
    slibreal_t x0, slibreal_t y0, slibreal_t z,
    slibreal_t px0, slibreal_t py0
) {
    Vector<3> detx = detector->GetPosition(), rcp;
    bool has_outer_wall = (wall_opacity!=WALL_OPACITY_SEMI_TRANSPARENT);

    for (unsigned int j = 0; j < ntoroidal; j++) {
        slibreal_t x = x0*cosphi[j] + y0*sinphi[j];
        slibreal_t y =-x0*sinphi[j] + y0*cosphi[j];

        slibreal_t px = px0*cosphi[j] + py0*sinphi[j];
        slibreal_t py =-px0*sinphi[j] + py0*cosphi[j];

        // Check if within FOV
        if (!IsWithinFieldOfView(x, y, z, rcp))
            continue;

        rp.UpdateXY(x, y, px, py, rcp);

        model->HandleParticle(&rp, otype, sinphi[j], cosphi[j]);

        if (!model->IsNonZero()) continue;
        else if (wall_opacity == WALL_OPACITY_TRANSPARENT ||
            !magfield->IntersectsDomain3D(
            x, y, z, detx[0], detx[1], detx[2], has_outer_wall)
        ) {
            RegisterOutput(&rp);
        }
    }
}
