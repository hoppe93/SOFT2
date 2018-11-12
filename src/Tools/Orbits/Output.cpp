/**
 * Implementation of output routines for the
 * 'Orbits' tool.
 */

#include <softlib/SFile.h>
#include "Tools/Orbits.h"

/**
 * Output.
 */
void Orbits::Output() {
    SFile *sf = SFile::Create(this->output, SFILE_MODE_WRITE);

    sf->WriteArray("t", this->tau, this->norbits, this->ntau);
    sf->WriteArray("x", this->x, this->norbits, 3*this->ntau);
    sf->WriteArray("p", this->p, this->norbits, 3*this->ntau);
    sf->WriteArray("solution", this->solution, this->norbits, 6*this->ntau);

    if (this->computeJacobian) {
        sf->WriteArray("Jdtdrho", this->Jdtdrho, this->norbits, this->ntau);
    }

    sf->Close();
}
