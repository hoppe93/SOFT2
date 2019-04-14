/**
 * Implementation of the 'Green' output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Output/Green.h"

using namespace __Radiation;
using namespace std;

/**
 * Constructor.
 */
Green::Green(Detector *d, MagneticField2D *m, ParticleGenerator *pgen) : RadiationOutput(d, m, pgen) {
    this->nr = pgen->GetNr();
    this->n1 = pgen->GetN1();
    this->n2 = pgen->GetN2();

    this->rgrid = pgen->GetRGrid();
    this->p1grid = pgen->GetP1Grid();
    this->p2grid = pgen->GetP2Grid();
    this->p1type = pgen->GetP1Type();
    this->p2type = pgen->GetP2Type();

    this->mpi_distribute_mode = pgen->GetDistributeMode();
    this->end_r = pgen->GetEndR();
    this->end_1 = pgen->GetEnd1();
    this->end_2 = pgen->GetEnd2();

    this->start_r = pgen->GetStartR();
    this->start_1 = pgen->GetStart1();
    this->start_2 = pgen->GetStart2();

    this->nw = d->GetNWavelengths();
}

/**
 * Destructor.
 */
Green::~Green() {
    // NOTE: These are NOT pointing to the arrays
    // the ParticleGenerator object after 'PrepareAllocateGreen()'
    // has been called, and can be safely deleted.
    delete [] this->rgrid;
    delete [] this->p1grid;
    delete [] this->p2grid;
}

