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

    // Write classification
    int32_t *cl = new int32_t[this->norbits];
    for (unsigned int i = 0; i < this->norbits; i++)
        cl[i] = (int32_t)this->classification[i];

    sf->WriteInt32List("classification", cl, this->norbits);
    delete [] cl;

    /*slibreal_t *cl = new slibreal_t[this->norbits];
    for (unsigned int i = 0; i < this->norbits; i++)
        cl[i] = (slibreal_t)this->classification[i];

    sf->WriteList("classification", cl, this->norbits);*/

    if (this->computeJacobian) {
        sf->WriteArray("Jdtdrho", this->Jdtdrho, this->norbits, this->ntau);
    }

    sf->Close();
}
