/**
 * Implementation of the sovvolume output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Output/SoVVolume.h"

using namespace __Radiation;

/**
 * Constructor.
 */
SoVVolume::SoVVolume(Detector *d, MagneticField2D *m, ParticleGenerator *pgen)
    : RadiationOutput(d, m) {

    this->np1    = pgen->GetN1();
    this->np2    = pgen->GetN2();
    this->p1grid = pgen->GetP1Grid();
    this->p2grid = pgen->GetP2Grid();
    this->p1type = pgen->GetP1Type();
    this->p2type = pgen->GetP2Type();
}

/**
 * Destructor.
 */
SoVVolume::~SoVVolume() { }

