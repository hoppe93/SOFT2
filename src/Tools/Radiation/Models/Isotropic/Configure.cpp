/**
 * Configuration of the 'Isotropic' radiation model.
 */

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Tools/Radiation/Models/Isotropic.h"

using namespace std;
using namespace __Radiation;

/**
 * Configure the isotropic radiation model.
 *
 * conf: ConfigBlock for the model.
 * root: Root ConfigBlock. Provided to allow access
 *       to configuration of possible sub-modules.
 *       Not currently used in the 'isotropic' model.
 */
void Isotropic::Configure(struct global_settings *__UNUSED__(globset), ConfigBlock *conf, ConfigBlock *__UNUSED__(root)) {
    // value
    if (conf->HasSetting("value")) {
        Setting *s = conf->GetSetting("value");
        if (!s->IsScalar())
            throw new IsotropicException("Invalid value assigned to parameter 'value'. Expected real number.");

        this->value = s->GetScalar();
    } else this->value = 1.0;
}

/**
 * Returns a string describing this model.
 */
const string Isotropic::GetDescription() const {
#ifdef COLOR_TERMINAL
    return Model::ISOTROPIC_COLOR+"Isotropic\x1B[0m";
#else
    return "Isotropic";
#endif
}

/**
 * Set the gyro-averaged amount of radiation emitted
 * by each particle per unit time and solid angle.
 */
void Isotropic::SetRadiationValue(slibreal_t P0) {
    this->value = P0;
}

