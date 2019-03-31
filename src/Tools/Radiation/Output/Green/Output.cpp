/**
 * Radiation :: Output :: Green
 *
 * Generate an output Green's function.
 */

#include <omp.h>
#include <softlib/config.h>
#include "Tools/Radiation/Output/Green.h"

using namespace __Radiation;

slibreal_t *Green::global_function = nullptr;

/**
 * Finish gathering output and merge output
 * stored on individual threads to the global
 * Green's function.
 */
void Green::Finish() {
    #pragma omp critical (Green_Finish)
    {
        if (global_function == nullptr) {
            global_function = this->function;
        } else {
            for (unsigned int i = 0; i < this->fsize; i++)
                global_function[i] += this->function[i];

            delete [] this->function;
        }
    }
}

/**
 * Generate the output file.
 * Called on the root thread only.
 */
void Green::Generate() {
    SFile *sf = SFile::Create(this->output, SFILE_MODE_WRITE);

    sf->WriteList("func", global_function, fsize);
    sf->WriteList("r", this->rgrid, this->nr);
    sf->WriteList("param1", this->p1grid, this->n1);
    sf->WriteList("param2", this->p2grid, this->n2);
    sf->WriteString("param1name", Particle::GetCoordinateName(this->p1type));
    sf->WriteString("param2name", Particle::GetCoordinateName(this->p2type));
    sf->WriteString("type", this->format);
    sf->WriteList("wavelengths", this->detector->GetWavelengths(), this->nw);

    if (this->hasI)
        sf->WriteScalar("rowpixels", this->nrowpixels);
    if (this->hasJ)
        sf->WriteScalar("colpixels", this->ncolpixels);

    if (this->storeStokesParameters)
        sf->WriteScalar("stokesparams", 1.0);
	else
		sf->WriteScalar("stokesparams", 0.0);

	this->WriteCommonQuantities(sf);

    sf->Close();

    SOFT::PrintInfo("Wrote Green's function to '%s'.", this->output.c_str());

    delete [] global_function;
}

