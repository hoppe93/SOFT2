/**
 * Initialization of SOFT
 *
 * This module reads the configuration and creates a new SOFT object,
 * populating it with the necessary objects.
 */

#include <iostream>
using namespace std;

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>

#include "config.h"
#include "Init/Init.h"
#include "Init/InitGlobal.h"
#include "Init/InitMagneticField.h"
#include "Init/InitPusher.h"
#include "PhaseSpace/ParticleGenerator.h"
#include "SOFT.h"
#include "SOFTException.h"

const string soft_init_defaults=
"distribution_function=__unit_distribution_function__;\n"
"include_drifts=no;\n"
"magnetic_field=__default__;\n"
"num_threads=__default__;\n"
"particle_generator=__default__;\n"
"particle_pusher=__default__;\n"
"@DistributionFunction __unit_distribution_function__ (unit) {}"
;

/**
 * Create a new SOFT object and initialize it
 * according to the given configuration.
 *
 * conf: Configuration object specifying how to
 *       initialize the SOFT object and its children.
 */
SOFT *InitSOFT(Configuration *input) {
    Configuration *conf;
    ConfigBlock *cfb, root;
    SOFT *soft;
    string defaults;

    soft = new SOFT();

    defaults = InitSOFTDefaults();

    // Override default settings
    conf = new Configuration();
    InitConfig(conf);
    conf->FromString(defaults, "<SOFT defaults>");
    conf->Merge(*input, true);

    soft->configuration = conf;

    // GLOBAL SETTINGS
    root = conf->GetRootBlock();
    soft->SetGlobalSettings(InitGlobalSettings(root));

    // MAGNETIC FIELD
    cfb = root.GetConfigBlock(CONFBLOCK_MAGNETICFIELD, soft->GetGlobalSettings()->magnetic_field);
    if (cfb == nullptr)
        throw SOFTException("Magnetic field '%s' not defined.", soft->GetGlobalSettings()->magnetic_field.c_str());
    soft->magfield = InitMagneticField(cfb);

    // DISTRIBUTION FUNCTION
    cfb = root.GetConfigBlock(CONFBLOCK_DISTRIBUTION, soft->GetGlobalSettings()->distribution);
    if (cfb == nullptr)
        throw SOFTException("Distribution function '%s' not defined.", soft->GetGlobalSettings()->distribution.c_str());
    soft->distribution = InitDistributionFunction(soft->magfield, cfb, &root);
    
    // PARTICLE GENERATOR
    cfb = root.GetConfigBlock(CONFBLOCK_PARTICLEGENERATOR, soft->GetGlobalSettings()->particle_generator);
    if (cfb == nullptr)
        throw SOFTException("Particle generator '%s' not defined.", soft->GetGlobalSettings()->particle_generator.c_str());
    soft->partgen = new ParticleGenerator(soft->magfield, cfb, soft->GetGlobalSettings());

    // Limit number of threads if necessary
    if (soft->GetGlobalSettings()->num_threads > soft->partgen->Size())
        soft->GetGlobalSettings()->num_threads = soft->partgen->Size();
    
    return soft;
}

/**
 * Create a default configuration.
 */
string InitSOFTDefaults() {
    string defaults = soft_init_defaults;

    defaults += InitPusherDefaults();
    return defaults;
}

/*********************************
 * HELPERS                       *
 *********************************/
/**
 * Tries to read the parameter named 'name' from
 * the ConfigBlock 'conf' as a numeric scalar value.
 * If it fails (i.e. if the parameter does not exist
 * or if the assigned value is not a scalar), a SOFTException
 * is thrown which contains information about the
 * parameter as well as which parent block (named 'parent')
 * that contained the error.
 *
 * conf:   ConfigBlock supposed to contain the parameter.
 * name:   Name of parameter to read.
 * parent: Name of parent block to read parameter from.
 */
slibreal_t init_get_scalar(ConfigBlock *conf, const string& name, const string& parent) {
    Setting *s;
    if (!conf->HasSetting(name))
        throw SOFTException("%s: Required parameter '%s' not defined.", parent.c_str(), name.c_str());
    
    s = conf->GetSetting(name);
    if (!s->IsScalar())
        throw SOFTException("%s: %s: Invalid value assigned to parameter. Expected real scalar value.", parent.c_str(), name.c_str());

    return s->GetScalar();
}

/**
 * Tries to read the parameter 'name' from
 * the ConfigBlock 'conf' as an integer. If it
 * fails (i.e. if the parameter does not exist
 * or if the assigned value is not an integer),
 * a SOFTException is thrown which contains
 * information about the parameter as well as which
 * parent block (named 'parent') that contained
 * the error.
 *
 * conf:   ConfigBlock supposed to contain parameter.
 * name:   Name of parameter to read.
 * parent: Name of the parent block to read parameter from.
 */
uint32_t init_get_uint32(ConfigBlock *conf, const string& name, const string& parent) {
    Setting *s;
    if (!conf->HasSetting(name))
        throw SOFTException("%s: Required parameter '%s' not defined.", parent.c_str(), name.c_str());

    s = conf->GetSetting(name);
    if (!s->IsUnsignedInteger32())
        throw SOFTException("%s: %s: Invalid value assigned to parameter. Expected non-negative 32-bit integer.", parent.c_str(), name.c_str());

    return s->GetUnsignedInteger32();
}

/**
 * Tries to read the parameter named 'name' from
 * the ConfigBlock 'conf' as a single string value.
 * If it fails (i.e. if the parameter does not exist
 * or if there are multiple assigned values), a SOFTException
 * is thrown which contains information about the
 * parameter as well as which parent block (named 'parent')
 * that contained the error.
 *
 * conf:   ConfigBlock supposed to contain the parameter.
 * name:   Name of parameter to read.
 * parent: Name of parent block to read parameter from.
 */
string& init_get_string(ConfigBlock *conf, const string& name, const string& parent) {
    Setting *s;
    if (!conf->HasSetting(name))
        throw SOFTException("%s: Required parameter '%s' not defined.", parent.c_str(), name.c_str());

    s = conf->GetSetting(name);
    if (s->GetNumberOfValues() != 1)
        throw SOFTException("%s: %s: Invalid value assigned to parameter. Expected single string value.", parent.c_str(), name.c_str());

    return s->GetString();
}

