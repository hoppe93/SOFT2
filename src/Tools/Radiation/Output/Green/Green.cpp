/**
 * Implementation of the 'Green' output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Output/Green.h"

using namespace __Radiation;

/**
 * Constructor.
 */
Green::Green(Detector *d, MagneticField2D *m, ParticleGenerator *pgen) : RadiationOutput(d, m) {
    this->nr = pgen->GetNr();
    this->n1 = pgen->GetN1();
    this->n2 = pgen->GetN2();
    this->rgrid = pgen->GetRGrid();
    this->p1grid = pgen->GetP1Grid();
    this->p2grid = pgen->GetP2Grid();
    this->p1type = pgen->GetP1Type();
    this->p2type = pgen->GetP2Type();

    this->nw = d->GetNWavelengths();
}

/**
 * Destructor.
 */
Green::~Green() { }

