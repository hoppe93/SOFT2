/**
 * Radiation :: Output :: Spectrum
 *
 * Configuration of the 'Spectrum' radiation output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Tools/Radiation/Output/Spectrum.h"

using namespace __Radiation;

/**
 * Configure the 'Spectrum' radiation output module.
 * 
 * conf: Configuration block for the module.
 * root: Root configuration block, providing access
 *       to the ConfigBlock's of possible sub-modules (unused).
 */
void Spectrum::Configure(ConfigBlock *conf, ConfigBlock *__UNUSED__(root)) {
    this->SetName(conf->GetName());

    // output
    if (!conf->HasSetting("output"))
        throw SpectrumException("No spectrum output name has been specifed. Set the option 'output'.");
    else
        this->output = (*conf)["output"];
}

/**
 * Initialize this output after configuration,
 * but before being used.
 */
void Spectrum::Initialize() { }

