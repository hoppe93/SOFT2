/**
 * Configuration of an 'Orbits' object.
 */

#include <string>
#include <softlib/Configuration.h>
#include "Tools/Orbits.h"


const std::string Orbits::DEFAULT_QUANTITIES[] = {
    OutputModule::PARAM1,
    OutputModule::PARAM1NAME,
    OutputModule::PARAM2,
    OutputModule::PARAM2NAME,
    OutputModule::R,
	OutputModule::RO_WALL
};
template<typename T, unsigned int sz>
unsigned int __def_size(T(&)[sz]) { return sz; }
const unsigned int Orbits::NDEFAULT_QUANTITIES = __def_size(Orbits::DEFAULT_QUANTITIES);

/**
 * Configure the 'Orbits' object.
 *
 * conf: Configuration object containing settings.
 */
void Orbits::Configure(
    struct global_settings *__UNUSED__(globset),
    ConfigBlock *conf, ConfigBlock *__UNUSED__(root)
) {
    this->Tool::SetName(conf->GetName());
    this->OutputModule::SetName(conf->GetName());

	// common quantities
	if (conf->HasSetting("common"))
		this->ConfigureCommonQuantities(DEFAULT_QUANTITIES, NDEFAULT_QUANTITIES, conf->GetSetting("common")->GetTextVector());
	else
		this->ConfigureCommonQuantities(DEFAULT_QUANTITIES, NDEFAULT_QUANTITIES);

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
