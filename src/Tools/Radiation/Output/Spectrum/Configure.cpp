/**
 * Radiation :: Output :: Spectrum
 *
 * Configuration of the 'Spectrum' radiation output module.
 */

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Tools/Radiation/Output/Spectrum.h"

using namespace __Radiation;
using namespace std;

const string Spectrum::DEFAULT_QUANTITIES[] = { "none" };
/*template<typename T, unsigned int sz>
unsigned int __def_size(T(&)[sz]) { return sz; }
const unsigned int Spectrum::NDEFAULT_QUANTITIES = __def_size(Spectrum::DEFAULT_QUANTITIES);*/
const unsigned int Spectrum::NDEFAULT_QUANTITIES = 0;

/**
 * Configure the 'Spectrum' radiation output module.
 * 
 * conf: Configuration block for the module.
 * root: Root configuration block, providing access
 *       to the ConfigBlock's of possible sub-modules (unused).
 */
void Spectrum::Configure(ConfigBlock *conf, ConfigBlock *__UNUSED__(root)) {
    this->SetName(conf->GetName());

	// common
	if (conf->HasSetting("common"))
		this->ConfigureCommonQuantities(DEFAULT_QUANTITIES, NDEFAULT_QUANTITIES, conf->GetSetting("common")->GetTextVector());
	else
		this->ConfigureCommonQuantities(DEFAULT_QUANTITIES, NDEFAULT_QUANTITIES);

    // output
    if (!conf->HasSetting("output"))
        throw SpectrumException("No spectrum output name has been specifed. Set the option 'output'.");
    else
        this->output = (*conf)["output"];

    // stokesparams
    if (conf->HasSetting("stokesparams")) {
        Setting *s = conf->GetSetting("stokesparams");

        if (!s->IsBool())
            throw SpectrumException("Illegal value assigned to 'stokesparams'. Expected boolean.");

        this->measurePolarization = s->GetBool();
    }
}

/**
 * Initialize this output after configuration,
 * but before being used.
 */
void Spectrum::Initialize() {
	AllocateSpectrum();
}

void Spectrum::AllocateSpectrum() {
    this->I = new slibreal_t[this->nwavelengths];

    for (unsigned int i = 0; i < this->nwavelengths; i++)
        this->I[i] = 0.0;

    if (this->MeasuresPolarization()) {
        this->Q = new slibreal_t[this->nwavelengths];
        this->U = new slibreal_t[this->nwavelengths];
        this->V = new slibreal_t[this->nwavelengths];

        for (unsigned int i = 0; i < this->nwavelengths; i++) {
            this->Q[i] = 0.0;
            this->U[i] = 0.0;
            this->V[i] = 0.0;
        }
    }
}

