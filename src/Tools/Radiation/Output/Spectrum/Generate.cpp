/**
 * Radiation :: Output :: Spectrum
 *
 * Generate an output spectrum file.
 */

#include <omp.h>
#include <softlib/config.h>
#include "Tools/Radiation/Output/Spectrum.h"

using namespace __Radiation;

slibreal_t
    *Spectrum::global_I = nullptr,
    *Spectrum::global_Q = nullptr,
    *Spectrum::global_U = nullptr,
    *Spectrum::global_V = nullptr;

/**
 * Finish gathering output and merge output
 * stored on individual threads to the global
 * spectrum.
 */
void Spectrum::Finish() {
    #pragma omp critical (Spectrum_Finish)
    {
        if (global_I == nullptr) {
            global_I = this->I;
            if (this->MeasuresPolarization()) {
                global_Q = this->Q;
                global_U = this->U;
                global_V = this->V;
            }
        } else {
            for (unsigned int i = 0; i < this->nwavelengths; i++)
                global_I[i] += this->I[i];

            delete [] this->I;

            if (this->MeasuresPolarization()) {
                for (unsigned int i = 0; i < this->nwavelengths; i++) {
                    global_Q[i] += this->Q[i];
                    global_U[i] += this->U[i];
                    global_V[i] += this->V[i];
                }

                delete [] this->Q;
                delete [] this->U;
                delete [] this->V;
            }
        }
    }
}

/**
 * Generate the output file.
 * Called on the root thread only.
 */
void Spectrum::Generate() {
    SFile *sf = SFile::Create(this->output, SFILE_MODE_WRITE);

    sf->WriteList("wavelengths", this->wavelengths, this->nwavelengths);
    sf->WriteList("I", global_I, this->nwavelengths);
    if (this->MeasuresPolarization()) {
        sf->WriteList("Q", global_Q, this->nwavelengths);
        sf->WriteList("U", global_U, this->nwavelengths);
        sf->WriteList("V", global_V, this->nwavelengths);
    }

	this->WriteCommonQuantities(sf);

    sf->Close();

    SOFT::PrintInfo("Wrote spectrum to '%s'.", this->output.c_str());

    if (this->MeasuresPolarization()) {
        delete [] global_V;
        delete [] global_U;
        delete [] global_Q;
    }

    delete [] global_I;
}

