/**
 * Configuration of the 'AngularDistribution' radiation model.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "SOFT.h"
#include "Tools/Radiation/Models/AngularDistribution.h"
#include "Tools/Radiation/Models/AngularDistribution/ADSynchrotronEmission.h"
#include "Tools/Radiation/Models/AngularDistribution/Quadrature2D/ADEval2D.h"
#include "Tools/Radiation/Models/AngularDistribution/Quadrature2D/ADSimpson2D.h"

using namespace std;
using namespace __Radiation;

const unsigned int AngularDistribution::NSAMPLES_DEFAULT=1;

/**
 * Configure the angular distribution model.
 *
 * conf: ConfigBlock for the model.
 * root: Root ConfigBlock. Provided to allow access
 *       to configuration of possible sub-modules.
 */
void AngularDistribution::Configure(struct global_settings *globset, ConfigBlock *conf, ConfigBlock *__UNUSED__(root)) {
    // emission  -- Name of emission model to use
    if (!conf->HasSetting("emission"))
        throw AngularDistributionException("No emission type specified.");
    ConfigureEmission(conf, (*conf)["emission"], globset);

    // nsamples  -- Number of sample points on detector
    if (conf->HasSetting("nsamples")) {
        Setting *s = conf->GetSetting("nsamples");

        if (!s->IsUnsignedInteger32() || s->GetUnsignedInteger32() == 0)
            throw AngularDistributionException("Invalid value assigned to 'nsamples'. Expected positive integer.");
        this->nsamples = s->GetUnsignedInteger32();
    } else
        this->nsamples = NSAMPLES_DEFAULT;

    // qrule2d  -- Quadruature rule to use for the 2D integral
    //             (must be done after 'emission')
    if (nsamples == 1) {
        if (conf->HasSetting("qrule2D"))
            SOFT::PrintInfo("Number of samples on detector surface is 1. Ignoring setting 'qrule2d'.");
        this->quadrature2d = new ADEval2D(this->emission, this->parent->detector);
    } else if (conf->HasSetting("qrule2d")) {
        if ((*conf)["qrule2d"] == "simpson")
            this->quadrature2d = new ADSimpson2D(this->emission, this->parent->detector, this->nsamples);
        else
            throw AngularDistributionException("Invalid quadrature rule specified.");
    } else
        this->quadrature2d = new ADSimpson2D(this->emission, this->parent->detector, this->nsamples);
}

/**
 * Configure the emission module.
 *
 * emname: Name of the emission module to use.
 */
void AngularDistribution::ConfigureEmission(ConfigBlock *conf, const string& emname, struct global_settings *globset) {
    if (emname == "bremsstrahlung") {
        throw AngularDistributionException("Bremsstrahlung has not been implemented in the angular distribution model yet.");
    } else if (emname == "synchrotron") {
        this->emission = ConfigureSynchrotronEmission(conf, globset);
    } else
        throw AngularDistributionException("Unrecognized emission model requested: '%s'.", emname.c_str());

    this->emissionName = emname;
}

/**
 * Configure synchrotron emission.
 */
ADSynchrotronEmission *AngularDistribution::ConfigureSynchrotronEmission(ConfigBlock *conf, struct global_settings *globset) {
    Setting *set;
    size_t qagsLimit = 100;
    slibreal_t qagsEpsRel = 1e-3;

    if (conf->HasSetting("qagslimit")) {
        set = conf->GetSetting("qagslimit");
        if (!set->IsUnsignedInteger64())
            throw AngularDistributionException("Invalid value assigned to 'qagstol'. Expected 64-bit unsigned integer.");

        qagsLimit = set->GetUnsignedInteger64();
    }

    if (conf->HasSetting("qagstol")) {
        set = conf->GetSetting("qagstol");
        if (!set->IsScalar())
            throw AngularDistributionException("Invalid value assigned to 'qagstol'. Expected real scalar.");

        qagsEpsRel = set->GetScalar();
    }

    return new ADSynchrotronEmission(
        this->parent->detector, globset,
        qagsLimit, qagsEpsRel
    );
}

/**
 * Returns a string describing this model.
 */
const string AngularDistribution::GetDescription() const {
#ifdef COLOR_TERMINAL
    string model = Model::ANGDIST_COLOR+"Angular distribution\e[0m  with ";
    if (this->emissionName == "bremsstrahlung")
        model += Model::BREMSSTRAHLUNG_COLOR+"bremsstrahlung\e[0m";
    else if (this->emissionName == "synchrotron")
        model += Model::SYNCHROTRON_COLOR+"synchrotron\e[0m";
    else
        model += Model::UNIT_EMISSION_COLOR+this->emissionName+"\e[0m";

    return model;
#else
    return "Angular distribution  with "+this->emissionName;
#endif
}
