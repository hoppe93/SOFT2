/**
 * Configuration of the 'Cone' radiation model.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "SOFT.h"
#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Cone/ConeBremsstrahlungEmission.h"
#include "Tools/Radiation/Models/Cone/ConeSynchrotronEmission.h"
#include "Tools/Radiation/Models/Cone/ConeUnitEmission.h"
#include "Tools/Radiation/Models/Cone/Projection/ConeProjection.h"
#include "Tools/Radiation/Models/Cone/Projection/Original.h"
#include "Tools/Radiation/Models/Cone/Projection/Reverse.h"

using namespace std;
using namespace __Radiation;

/**
 * Configure the cone model.
 *
 * conf: ConfigBlock for the model.
 * root: Root ConfigBlock. Provided to allow access
 *       to configuration of possible sub-modules.
 */
void Cone::Configure(struct global_settings *__UNUSED__(globset), ConfigBlock *conf, ConfigBlock *__UNUSED__(root)) {
    // edgecheck
    if (conf->HasSetting("edgecheck")) {
        Setting *s = conf->GetSetting("edgecheck");
        if (!s->IsBool())
            throw ConeException("Invalid value assigned to parameter 'edgecheck'. Expected boolean.");

        this->edgeCheck = s->GetBool();
    } else
        this->edgeCheck = false;

    // zeff (must be set before 'emission' model
    slibreal_t Zeff = 1;
    if (conf->HasSetting("zeff")) {
        Setting *s = conf->GetSetting("zeff");
        if (!s->IsScalar())
            throw ConeException("Invalid value assigned to paramter 'zeff'. Expected real scalar.");

        Zeff = s->GetScalar();
    }

    // emission
    if (!conf->HasSetting("emission"))
        throw ConeException("No emission type specified.");
    ConfigureEmission((*conf)["emission"], Zeff);

    // projection
    if (conf->HasSetting("projection")) {
        if ((*conf)["projection"] == "original") {
            this->projection = new ConeProjectionOriginal(this->parent->detector);
            SOFT::PrintWarning(SOFT::WARNING_TRMC_BUGGY_CONE_MODEL, "A buggy cone model is being used. Please, be careful.");
        } else if ((*conf)["projection"] == "reverse")
            this->projection = new ConeProjectionReverse(this->parent->detector);
        else
            throw ConeException("Invalid value assigned to parameter 'projection'.");
    } else
        this->projection = new ConeProjectionReverse(this->parent->detector);
}

/**
 * Configure the emission module.
 *
 * emname: Name of the emission module to use.
 */
void Cone::ConfigureEmission(const string& emname, const slibreal_t zeff) {
    if (emname == "bremsstrahlung") {
        if (this->parent->MeasuresPolarization())
            throw ConeBremsstrahlungException("The bremsstrahlung radiation model does not support polarization measurements.");

        this->emission = new ConeBremsstrahlungEmission(this->parent->detector, zeff);
    } else if (emname == "synchrotron") {
        this->emission = new ConeSynchrotronEmission(this->parent->detector);
    } else if (emname == "unit") {
        if (this->parent->MeasuresPolarization())
            throw ConeUnitException("The 'Unit' emission model does not support polarization measurements.");

        this->emission = new ConeUnitEmission(this->parent->detector);
    } else
        throw ConeException("Unrecognized emission model requested: '%s'.", emname.c_str());

    this->emissionName = emname;
}

/**
 * Returns a string describing this model.
 */
const string Cone::GetDescription() const {
#ifdef COLOR_TERMINAL
    string model = Model::CONE_COLOR+"Cone\x1B[0m  with ";
    if (this->emissionName == "bremsstrahlung")
        model += Model::BREMSSTRAHLUNG_COLOR+"bremsstrahlung\x1B[0m";
    else if (this->emissionName == "synchrotron")
        model += Model::SYNCHROTRON_COLOR+"synchrotron\x1B[0m";
    else
        model += Model::UNIT_EMISSION_COLOR+this->emissionName+"\x1B[0m";

    return model;
#else
    return "Cone  with "+this->emissionName;
#endif
}

