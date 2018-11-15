/**
 * Initialization of global settings.
 */

#include <omp.h>
#include <vector>

#include <softlib/config.h>
#include <softlib/Configuration.h>

#include "Init/Init.h"
#include "SOFT.h"
#include "SOFTException.h"

using namespace std;

/**
 * Initialize a global settings object given
 * the root configuration block.
 *
 * global: Configuration block containing the
 *         "global" settings.
 */
struct global_settings *InitGlobalSettings(ConfigBlock& global) {
    struct global_settings *globset = new struct global_settings;

    // INCLUDE DRIFTS
    if (!global.GetSetting("include_drifts")->IsBool())
        throw SOFTException("include_drifts: Invalid value assigned to parameter. Expected boolean value.");
    else
        globset->include_drifts = global.GetSetting("include_drifts")->GetBool();

    // MAGNETIC FIELD
    if (global.GetSetting("magnetic_field")->GetNumberOfValues() != 1)
        throw SOFTException("magnetic_field: Invalid value assigned to parameter. Expected single name of magnetic field to use.");
    globset->magnetic_field = global.GetSetting("magnetic_field")->GetString();
    if (globset->magnetic_field == "__default__") {
        vector<ConfigBlock> cfb = global.GetAllSubBlocksOfType(CONFBLOCK_MAGNETICFIELD);
        if (cfb.size() == 0)
            throw SOFTException("No magnetic field defined.");
        else if (cfb.size() > 1)
            throw SOFTException("magnetic_field: More than one magnetic field defined, but which magnetic field to use was not specified.");
        else
            globset->magnetic_field = cfb[0].GetName();
    }

    // DISTRIBUTION FUNCTION
    if (global.GetSetting("distribution_function")->GetNumberOfValues() != 1)
        throw SOFTException("distribution_function: Invalid value assigned to parameter. Expected single name of distribution function to use.");
    globset->distribution = global.GetSetting("distribution_function")->GetString();

    // NUMBER OF THREADS
    if (!global.GetSetting("num_threads")->IsUnsignedInteger32()) {
        if (global.GetSetting("num_threads")->GetNumberOfValues() != 1 ||
            global.GetSetting("num_threads")->GetString() != "__default__")
            throw SOFTException("num_threads: Invalid value assigned to parameter. Expected single integer.");
        else
            globset->num_threads = omp_get_max_threads();
    } else {
        globset->num_threads = global.GetSetting("num_threads")->GetUnsignedInteger32();

        if (globset->num_threads == 0)
            throw SOFTException("num_threads: Invalid number of worker threads specified: %u\n", globset->num_threads);
    }

    // PARTICLE GENERATOR
    if (global.GetSetting("particle_generator")->GetNumberOfValues() != 1)
        throw SOFTException("particle_generator: Invalid value assigned to parameter. Expected single name of generator to use.");
    globset->particle_generator = global.GetSetting("particle_generator")->GetString();
    if (globset->particle_generator == "__default__") {
        vector<ConfigBlock> cfb = global.GetAllSubBlocksOfType(CONFBLOCK_PARTICLEGENERATOR);
        if (cfb.size() == 0)
            throw SOFTException("No particle generator defined.");
        else if (cfb.size() > 1)
            throw SOFTException("particle_generator: More than one particle generator defined, but which particle generator to use was not specified.");
        else
            globset->particle_generator = cfb[0].GetName();
    }

    // PARTICLE PUSHER
    if (global.GetSetting("particle_pusher")->GetNumberOfValues() != 1)
        throw SOFTException("particle_pusher: Invalid value assigned to parameter. Expected single name of generator to use.");
    globset->particle_pusher = global.GetSetting("particle_pusher")->GetString();
    if (globset->particle_pusher == "__default__") {
        vector<ConfigBlock> cfb = global.GetAllSubBlocksOfType(CONFBLOCK_PARTICLEPUSHER);
        if (cfb.size() == 0)
            throw SOFTException("No particle pusher defined.");
        else if (cfb.size() > 1)
            throw SOFTException("particle_pusher: More than one particle pusher defined, but which particle pusher to use was not specified.");
        else
            globset->particle_pusher = cfb[0].GetName();
    }

    // TOOLS
    if (!global.HasSetting("tools"))
        throw new SOFTException("No tools have been selected for the run.");
    globset->tools = global.GetSetting("tools")->GetTextVector();

    // Handle untouched settings
    if (global.HasUntouchedSettings()) {
        vector<Setting*> unt = global.GetUntouchedSettings();
        string err = "The following global settings were unrecognized: ";

        for (unsigned int i = 0; i < unt.size(); i++) {
            err += "'" + unt[i]->GetName() + "'";
            if (i+1 < unt.size())
                err += ", ";
        }

        err += ".";
        throw SOFTException("%s", err.c_str());
    }

    return globset;
}

