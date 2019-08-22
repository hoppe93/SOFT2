/**
 * Handle a single point of phase-space with the cone model.
 */

#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/Cone.h"

using namespace __Radiation;

/**
 * Constructor.
 */
Cone::Cone(Radiation *rad) : Model(rad) {
    unsigned int n = rad->detector->GetNWavelengths(), i;

    if (n > 0) {
        this->I = new slibreal_t[n];
        this->Q = new slibreal_t[n];
        this->U = new slibreal_t[n];
        this->V = new slibreal_t[n];

        for (i = 0; i < n; i++)
            I[i] = Q[i] = U[i] = V[i] = 0.0;
    } else
        this->I = this->Q = this->U = this->V = nullptr;
}
Cone::~Cone() {
    delete emission;
    delete projection;

    if (this->I != nullptr)
        delete [] this->I;
    if (this->Q != nullptr)
        delete [] this->Q;
    if (this->U != nullptr)
        delete [] this->U;
    if (this->V != nullptr)
        delete [] this->V;
}

/**
 * Calculates the amount of emitted radiation in this
 * timestep.
 *
 * rp: RadiationParticle defining the instantaneous
 *     particle properties.
 */
void Cone::InitializeTimestep(RadiationParticle*) {
    // Since the radiation polarization state depends on
    // a particle-frame coordinate system, we can't precompute
    // it (even though it is independent of the detector, it
    // is fixed in space and would have to be rotated, something
    // which is not easy/useful to implement)
    //if (!this->parent->MeasuresPolarization())
    //    this->emission->HandleParticle(rp, false);
}

/**
 * Handles the given particle radiation state and
 * calculates the fraction of the guiding-center cone
 * that overlaps the detector surface.
 *
 * rp: RadiationParticle defining the instantaneous
 *     particle properties.
 */
void Cone::HandleParticle(RadiationParticle *rp, orbit_type_t otype, const slibreal_t sinphi, const slibreal_t cosphi) {
    if (otype == ORBIT_TYPE_GUIDING_CENTER)
        ComputeOverlappingRadiationGC(rp, sinphi, cosphi);
    else if (otype == ORBIT_TYPE_PARTICLE)
        ComputeOverlappingRadiationParticle(rp); 
    else
        throw ConeException("Non-existent orbit type: %d", otype);
}

void Cone::ComputeOverlappingRadiationGC(RadiationParticle *rp, const slibreal_t sinphi, const slibreal_t cosphi){
    if (this->edgeCheck && !EdgeCheck(rp, sinphi, cosphi))
        overlapFraction = 0;
    else
        overlapFraction = projection->ComputeOverlappingFraction(rp);

    if (overlapFraction == 0) {
        this->nonzero = false;
    } else {
        this->nonzero = true;

        //if (this->parent->MeasuresPolarization())
            this->emission->HandleParticle(rp, this->parent->MeasuresPolarization());
    }

    // Multiply quantities by fraction
    this->totEmission = this->emission->GetTotalEmission() * overlapFraction;
    unsigned int n = this->emission->GetNWavelengths(), i;
    slibreal_t *s = this->emission->GetSpectrum();
    for (i = 0; i < n; i++)
        this->I[i] = s[i] * overlapFraction;

    if (this->parent->MeasuresPolarization()) {
        this->totQ = this->emission->GetPowerQ();
        this->totU = this->emission->GetPowerU();
        this->totV = this->emission->GetPowerV();

        s = this->emission->GetStokesQ();
        for (i = 0; i < n; i++)
            this->Q[i] = s[i] * overlapFraction;
        s = this->emission->GetStokesU();
        for (i = 0; i < n; i++)
            this->U[i] = s[i] * overlapFraction;

        // V is identically zero in the cone model!
    }
}

//Kollar om strålningen hamnar på detektorplanet
void Cone::ComputeOverlappingRadiationParticle(RadiationParticle *rp){ 
    Vector<3> Rdp = rp->GetRCP(),
        vhat = rp->GetPHat(),
        ehat1 = this->parent->detector->GetEHat1(),
        ehat2 = this->parent->detector->GetEHat2(),
        n = this->parent->detector->GetDirection();
    printf("vx = %.16e, vy = %.16e, vz = %.16e \n", Rdp[0], Rdp[1], Rdp[2]);
    slibreal_t HalfAp = this->parent->detector->GetAperture()/2;

    slibreal_t Rdp1 = Rdp.Dot(ehat1), Rdp2 = Rdp.Dot(ehat2), Rdpn = Rdp.Dot(n),
        vhat1 = vhat.Dot(ehat1), vhat2 =  vhat.Dot(ehat2), vhatn = vhat.Dot(n),
        a = -Rdpn/vhatn;

    if (fabs(Rdp1+a*vhat1) <= HalfAp && fabs(Rdp2+a*vhat2) <= HalfAp){
        this->nonzero = true;
        this->emission->HandleParticle(rp, this->parent->MeasuresPolarization());
        this->totEmission = this->emission->GetTotalEmission();
        unsigned int n = this->emission->GetNWavelengths(), i;
        slibreal_t *s = this->emission->GetSpectrum();
        for (i = 0; i < n; i++)
            this->I[i] = s[i];
    }
    else{
        //this->totEmission = 0;
        this->nonzero = false;
    }
    
}


