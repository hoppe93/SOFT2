/**
 * Radiation :: Output :: Topview
 *
 * Configuration of the 'Topview' radiation output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Tools/Radiation/Output/Topview.h"

using namespace __Radiation;

/**
 * Allocate memory for the topview.
 */
void Topview::AllocateTopview() {
    this->ntotpixels = this->npixels*this->npixels;
    this->topview = new slibreal_t[this->ntotpixels];

    // Initialize
    for (unsigned int i = 0; i < this->ntotpixels; i++)
        this->topview[i] = 0.0;
}

/**
 * Configure the 'Topview' radiation output module.
 * 
 * conf: Configuration block for the module.
 * root: Root configuration block, providing access
 *       to the ConfigBlock's of possible sub-modules (unused).
 */
void Topview::Configure(ConfigBlock *conf, ConfigBlock *__UNUSED__(root)) {
    Setting *s;

    this->SetName(conf->GetName());

    // output
    if (!conf->HasSetting("output"))
        throw TopviewException("No topview output name has been specifed. Set the option 'output'.");
    else
        this->output = (*conf)["output"];

    // pixels
    if (!conf->HasSetting("pixels"))
        throw TopviewException("The number of pixels in the topview was not set.");
    else {
        s = conf->GetSetting("pixels");
        if (s->GetNumberOfValues() == 1)
            this->npixels = s->GetUnsignedInteger32();
        else
            throw TopviewException("Invalid specification of 'pixels'.");
    }

    AllocateTopview();
}

