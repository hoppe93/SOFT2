/**
 * Configuration of an 'Orbits' object.
 */

#include <softlib/Configuration.h>
#include "Tools/Orbits.h"

/**
 * Configure the 'Orbits' object.
 *
 * conf: Configuration object containing settings.
 */
void Orbits::Configure(
    struct global_settings *__UNUSED__(globset),
    ConfigBlock *conf, ConfigBlock *__UNUSED__(root)
) {
    this->SetName(conf->GetName());

    if (!conf->HasSetting("output"))
        throw OrbitsException("Required parameter 'output' not defined.");
    else
        this->output = (*conf)["output"];

    if (conf->HasSetting("computeJacobian")) {
        Setting *s = conf->GetSetting("computeJacobian");
        if (!s->IsBool())
            throw OrbitsException("Unexpected value for parameter 'computeJacobian'.");
        ToggleJacobianCalculation(s->GetBool());
    } else
        ToggleJacobianCalculation(false);
}

/**
 * Prepare the SOFT configuration object to
 * read settings for this tool.
 *
 * conf: Configuration object to initialize.
 */
void Orbits::PrepareConfiguration(Configuration *__UNUSED__(conf)) { }
