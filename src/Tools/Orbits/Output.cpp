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
    sf->WriteArray("ppar", this->ppar, this->norbits, this->ntau);
    sf->WriteArray("pperp", this->pperp, this->norbits, this->ntau);
    sf->WriteArray("B", this->B, this->norbits, 3*this->ntau);
    sf->WriteArray("bhat", this->bhat, this->norbits, 3*this->ntau);
    sf->WriteArray("p2", this->p2, this->norbits, this->ntau);
    sf->WriteArray("Babs", this->Babs, this->norbits, this->ntau);
    sf->WriteArray("gamma", this->gamma, this->norbits, this->ntau);

    // Write classification
    int32_t *cl = new int32_t[this->norbits];
    for (unsigned int i = 0; i < this->norbits; i++)
        cl[i] = (int32_t)this->classification[i];

    sf->WriteInt32List("classification", cl, this->norbits);
    delete [] cl;

    sf->WriteList("driftshift", this->drift_shift, this->norbits);

    if (this->computeJacobian) {
        sf->WriteArray("Jdtdrho", this->Jdtdrho, this->norbits, this->ntau);
    }

    this->WriteCommonQuantities(sf);

    sf->Close();
}
