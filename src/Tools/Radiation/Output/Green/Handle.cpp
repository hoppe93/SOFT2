/**
 * Radiation :: Output :: Green
 *
 * Handle output from a RadiationModel.
 */

#include "config.h"
#include "Tools/Radiation/Output/Green.h"
#include "Tools/Radiation/Output/Image.h"

using namespace __Radiation;

/**
 * Handle the result from the given radiation model.
 *
 * det: Object representing the detector.
 * m:   Radiation model doing the computation.
 * rp:  Object representing the particle emitting state.
 */
void Green::Handle(Detector *det, Model *m, RadiationParticle *rp) {
    size_t index = 0, wavindex = 0;

    GetIndex(det, rp, &index, &wavindex);

    if (storeStokesParameters)
        index *= 4;

    // Compute differential element
    //slibreal_t diffel = rp->GetRDphi() * rp->GetJdtdrho() * rp->GetJp();
    slibreal_t diffel;

    if (this->withJacobian)
        diffel = rp->GetRDphi() * rp->GetJdtdrho() * rp->GetJp();
    else
        //diffel = rp->GetRDphi() * rp->GetJp();
        diffel = 1.0;

    if (weighWithDistribution) diffel *= rp->GetF();
    if (!hasP1) diffel *= rp->GetDP1();
    if (!hasP2) diffel *= rp->GetDP2();
    if (hasR)   diffel /= rp->GetDrho();

    // Set value of Green's function
    if (hasW) {
        slibreal_t *I, *Q=nullptr, *U=nullptr, *V=nullptr;
        const size_t inc = (this->iw==this->nformat-1?0:this->factors[this->iw]);

        I = m->GetStokesI();
        if (storeStokesParameters) {
            Q = m->GetStokesQ();
            U = m->GetStokesU();
            V = m->GetStokesV();
        }

        // Functions for copying either Stokes parameters
        // or spectrum
        auto copyStokes = [this,&I,&Q,&U,&V,&diffel,&inc,&index]() {
            int i, j;
            for (i = 0; i < this->nw; i++, index+=inc) {
                j = index + 4*i;

                #pragma omp atomic update
                this->function[j+0] += I[i] * diffel;
                #pragma omp atomic update
                this->function[j+1] += Q[i] * diffel;
                #pragma omp atomic update
                this->function[j+2] += U[i] * diffel;
                #pragma omp atomic update
                this->function[j+3] += V[i] * diffel;
            }
        };
        auto copySpec = [this,&I,&diffel,&inc,&index]() {
            int i;
            for (i = 0; i < nw; i++, index+=inc) {
                #pragma omp atomic update
                this->function[index+i] += I[i] * diffel;
            }
        };

        /* If one or more momentum-space parameters are
         * NOT in the Green's function, then it's possible
         * for multiple threads to want to write to the
         * same element in the Green's function (because
         * each thread is assigned an individual point
         * of phase-space) and therefore we must block
         * other threads from writing at the same time.
         */
        if (!containsAllMomentumSpaceParameters) {
            if (storeStokesParameters)
                copyStokes();
            else copySpec();
        } else {
            if (storeStokesParameters)
                copyStokes();
            else copySpec();
        }
    } else {
        slibreal_t I, Q=0, U=0, V=0;

        I = m->GetPower();
        if (storeStokesParameters) {
            Q = m->GetPowerQ();
            U = m->GetPowerU();
            V = m->GetPowerV();
        }

        /* If one or more momentum-space parameters are
         * NOT in the Green's function, then it's possible
         * for multiple threads to want to write to the
         * same element in the Green's function (because
         * each thread is assigned an individual point
         * of phase-space) and therefore we must block
         * other threads from writing at the same time.
         */
        if (!containsAllMomentumSpaceParameters) {
            if (storeStokesParameters) {
                #pragma omp atomic update
                this->function[index+0] += I * diffel;
                #pragma omp atomic update
                this->function[index+1] += Q * diffel;
                #pragma omp atomic update
                this->function[index+2] += U * diffel;
                #pragma omp atomic update
                this->function[index+3] += V * diffel;
            } else {
                #pragma omp atomic update
                this->function[index] += I * diffel;
            }
        } else {
            if (storeStokesParameters) {
                this->function[index+0] += I * diffel;
                this->function[index+1] += Q * diffel;
                this->function[index+2] += U * diffel;
                this->function[index+3] += V * diffel;
            } else
                this->function[index] += I * diffel;
        }
    }
}

/**
 * Computes the Green's function index corresponding
 * to the given particle.
 *
 * det: Object representing the detector used.
 * rp:  Particle emitting state.
 */
void Green::GetIndex(
    Detector *det, RadiationParticle *rp,
    size_t *index, size_t *wavindex
) {
    int I, J;

    // Check if outside subimage
    if (hasI || hasJ) {
        Image::GetImagePixel(det, rp, this->nrowpixels, this->ncolpixels, I, J);

        if (I < suboffseti || I >= suboffseti+subnrowpixels)
            return;
        if (J < suboffsetj || J >= suboffsetj+subncolpixels)
            return;

        I -= suboffseti;
        J -= suboffsetj;
    }

    for (int i = 0; i < this->nformat; i++) {
        switch (format[i]) {
            case '1': *index += rp->GetIndexP1() * factors[i]; break;
            case '2': *index += rp->GetIndexP2() * factors[i]; break;
            case 'i': *index += I * factors[i]; break;
            case 'j': *index += J * factors[i]; break;
            case 'r': *index += rp->GetIndexR() * factors[i]; break;
            case 'w': *wavindex = i; break;
            default:
                throw GreenException("Unrecognized character in Green's function format: %c.", format[i]);
        }
    }
}

