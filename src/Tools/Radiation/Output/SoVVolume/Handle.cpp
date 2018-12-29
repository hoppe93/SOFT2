/**
 * Radiation :: Output :: SoVVolume
 * 
 * Handle output from a RadiationModel.
 */

#include <softlib/config.h>
#include "Tools/Radiation/Output/SoVVolume.h"

using namespace __Radiation;

/**
 * Handle output from the given model.
 *
 * det:   Detector observing the radiation.
 * model: Radiation model to handle output from.
 * rp:    Particle producing the radiation.
 */
void SoVVolume::Handle(Detector*, Model*, RadiationParticle *rp) {
    unsigned int
        i = rp->GetIndexP1(),
        j = rp->GetIndexP2();

    slibreal_t diffel = rp->GetRDphi() * rp->GetJdtdrho() * rp->GetJp();

    this->volumearray[i*this->np2 + j] += diffel;
}

