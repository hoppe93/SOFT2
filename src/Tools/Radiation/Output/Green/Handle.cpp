/**
 * Radiation :: Output :: Green
 *
 * Handle output from a RadiationModel.
 */

#include "config.h"
#include "Tools/Radiation/Output/Green.h"
#include "Tools/Radiation/Output/Image.h"

using namespace __Radiation;
using namespace std;

/**
 * Handle the result from the given radiation model.
 *
 * det: Object representing the detector.
 * m:   Radiation model doing the computation.
 * rp:  Object representing the particle emitting state.
 */
void Green::Handle(Detector *det, Model *m, RadiationParticle *rp) {
    size_t index = 0, wavindex = 0;

    // Get function index and make sure that the data
    // is located within the subimage
    if (!GetIndex(det, rp, &index, &wavindex))
        return;

    // Compute differential element
    //slibreal_t diffel = rp->GetRDphi() * rp->GetJdtdrho();
    slibreal_t diffel;

    if (this->withJacobian)
        diffel = rp->GetDphi() * rp->GetJdtdrho() * rp->GetDZeta();
    else
        diffel = 1.0;

    if (weighWithDistribution) diffel *= rp->GetF();
    if (!hasP1) diffel *= rp->GetDP1();
    if (!hasP2) diffel *= rp->GetDP2();
    if (hasR)   diffel /= rp->GetDrho();

    // Set value of Green's function
    if (hasW) {
        slibreal_t *I, *Q=nullptr, *U=nullptr, *V=nullptr;

        I = m->GetStokesI();
        if (storeStokesParameters) {
            Q = m->GetStokesQ();
            U = m->GetStokesU();
            V = m->GetStokesV();
        }

        // Functions for copying either Stokes parameters
        // or spectrum
        auto copyStokes = [this,&I,&Q,&U,&V,&diffel,&index,&wavindex]() {
            const size_t fsw = this->fsizeWithoutStokes;
            for (size_t i = 0; i < (unsigned)this->nw; i++) {
                #pragma omp atomic update
                this->function[fsw*0 + index + i*wavindex] += I[i] * diffel;
                #pragma omp atomic update
                this->function[fsw*1 + index + i*wavindex] += Q[i] * diffel;
                #pragma omp atomic update
                this->function[fsw*2 + index + i*wavindex] += U[i] * diffel;
                #pragma omp atomic update
                this->function[fsw*3 + index + i*wavindex] += V[i] * diffel;
            }
        };
        auto copySpec = [this,&I,&diffel,&index,&wavindex]() {
            for (int i = 0; i < nw; i++) {
                #pragma omp atomic update
                this->function[index+i*wavindex] += I[i] * diffel;
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
        if (!containsAllPhaseSpaceParameters) {
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
        if (!containsAllPhaseSpaceParameters) {
            if (storeStokesParameters) {
                #pragma omp atomic update
                this->function[this->fsizeWithoutStokes*0 + index] += I * diffel;
                #pragma omp atomic update
                this->function[this->fsizeWithoutStokes*1 + index] += Q * diffel;
                #pragma omp atomic update
                this->function[this->fsizeWithoutStokes*2 + index] += U * diffel;
                #pragma omp atomic update
                this->function[this->fsizeWithoutStokes*3 + index] += V * diffel;
            } else {
                #pragma omp atomic update
                this->function[index] += I * diffel;
            }
        } else {
            if (storeStokesParameters) {
                this->function[this->fsizeWithoutStokes*0 + index] += I * diffel;
                this->function[this->fsizeWithoutStokes*1 + index] += Q * diffel;
                this->function[this->fsizeWithoutStokes*2 + index] += U * diffel;
                this->function[this->fsizeWithoutStokes*3 + index] += V * diffel;
            } else
                this->function[index] += I * diffel;
        }
    }
}

/**
 * Computes the Green's function index corresponding
 * to the given particle. This function returns 'false'
 * if the given radiation particle lies outside the
 * subimage, and true otherwise.
 *
 * det: Object representing the detector used.
 * rp:  Particle emitting state.
 */
bool Green::GetIndex(
    Detector *det, RadiationParticle *rp,
    size_t *index, size_t *wavindex
) {
    int I, J;

    // Check if outside subimage
    if (hasI || hasJ) {
        Image::GetImagePixel(det, rp, this->nrowpixels, this->ncolpixels, I, J);

        if (I < suboffseti || I >= suboffseti+subnrowpixels)
            return false;
        if (J < suboffsetj || J >= suboffsetj+subncolpixels)
            return false;

        I -= suboffseti;
        J -= suboffsetj;
    }

    for (int i = 0; i < this->nformat; i++) {
        switch (format[i]) {
            case '1': *index += (rp->GetIndexP1() - this->start_1) * factors[i]; break;
            case '2': *index += (rp->GetIndexP2() - this->start_2) * factors[i]; break;
            case 'i': *index += I * factors[i]; break;
            case 'j': *index += J * factors[i]; break;
            case 'r': *index += (rp->GetIndexR() - this->start_r) * factors[i]; break;
            case 'w': *wavindex = this->factors[i]; break;
            default:
                throw GreenException("Unrecognized character in Green's function format: %c.", format[i]);
        }
    }

    return true;
}

