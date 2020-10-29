/**
 * Initialize tools and the ToolHandler.
 *
 * This should be done separately on each thread.
 */

#include <vector>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Init/Tools.h"
#include "SOFTException.h"
#include "SOFTLocal.h"
#include "Tools/Tool.h"
#include "Tools/ToolHandler.h"

#include "Tools/Integrator.h"
#include "Tools/Orbits.h"
#include "Tools/Radiation/Radiation.h"

using namespace std;

/***********************
 * TOOL INITIALIZATION *
 ***********************/
/**
 * HOW TO ADD A NEW TOOL
 *
 * 1. Add an entry to the 'toolinit_tools' array below.
 *    Each entry has four properties:
 *      [0]: Tool name (configuration block type name)
 *      [1]: Config block ID (set automatically, so can be anything)
 *      [2]: A pointer to the function 'InitTool' with an
 *           explicit type specification for the Tool class
 *           you're adding. For the 'Orbits' tool, this would be
 *
 *             InitTool<Orbits>
 *
 *      [3]: A pointer to the function 'PrepareToolConfiguration'
 *           with an explicit type specification for the Tool class
 *           you're adding. For the 'Orbits tool, this would be
 *
 *             PrepareToolConfiguration<Orbits>
 *
 * 2. Increase the 'INIT_TOOLS_NTOOLS' variable by one.
 */

template<class T>
Tool *InitTool(struct global_settings*, ConfigBlock*, ConfigBlock*, ParticleGenerator*, ParticlePusher*, MagneticField2D*);
template<class T>
void PrepareToolConfiguration(Configuration*);

/**
 * toolinit_tools array.
 * Add an entry for the new tool here. And
 * don't forget to increment 'INIT_TOOLS_NTOOLS' by one!
 */
const int INIT_TOOLS_NTOOLS = 3;
struct soft_toolspec toolinit_tools[INIT_TOOLS_NTOOLS] = {
    {"@Integrator", 0, InitTool<__SOFT::Integrator>, PrepareToolConfiguration<__SOFT::Integrator>},
    {"@Orbits",     0, InitTool<Orbits>,    PrepareToolConfiguration<Orbits>},
    {"@Radiation",  0, InitTool<__Radiation::Radiation>, PrepareToolConfiguration<__Radiation::Radiation>}
};

/**
 * Initialize a general Tool.
 *
 * globset: SOFT global settings.
 * conf:    Configuration block specifying the tool settings.
 * root:    Root configuration block.
 * pgen:    Particle generator used in the run.
 * pusher:  Particle pusher to work with.
 */
template<class T>
Tool *InitTool(
    struct global_settings *globset,
    ConfigBlock *conf, ConfigBlock *root,
    ParticleGenerator *pgen, ParticlePusher *pusher,
    MagneticField2D *mf
) {
    T *t = new T(mf, pgen, pusher);
    t->Configure(globset, conf, root);

    return t;
}

/**
 * Prepare for configuration of a Tool.
 *
 * conf: Configuration object to set up.
 */
template<class T>
void PrepareToolConfiguration(Configuration *conf) {
    T::PrepareConfiguration(conf);
}

/*******************************
 * TOOL HANDLER INITIALIZATION *
 *******************************/
/**
 * Initialize the ToolHandler and all used tools.
 *
 * globset: SOFT global settings.
 * input:   Configuration object to load settings from.
 * partgen: Particle generator object to use in run.
 */
ToolHandler *SOFTLocal::InitTools(
    struct global_settings *globset, Configuration *input,
    ParticleGenerator *partgen, ParticlePusher *pusher,
    MagneticField2D *magfield
) {
    unsigned int i, j;
    ToolHandler *th = new ToolHandler();
    ConfigBlock *root = input->GetRootBlock();

    vector<string> tools = globset->tools;
    for (i = 0; i < tools.size(); i++) {
        for (j = 0; j < INIT_TOOLS_NTOOLS; j++) {
            if (root->HasSubBlock(toolinit_tools[j].id, tools[i])) {
                ConfigBlock *cb = root->GetConfigBlock(toolinit_tools[j].id, tools[i]);
                Tool *t = toolinit_tools[j].init(globset, cb, root, partgen, pusher, magfield);

                th->AddTool(t);
                break;
            }
        }

        if (j == INIT_TOOLS_NTOOLS)
            throw SOFTException("Tool '%s' selected for run but not defined.", tools[i].c_str());
    }

    return th;
}

/**
 * Prepare the configuration handler to parse
 * configurations for tools.
 * 
 * conf: Configuration object to prepare.
 */
void InitToolConfig(Configuration *conf) {
    for (int i = 0; i < INIT_TOOLS_NTOOLS; i++) {
        toolinit_tools[i].id = conf->RegisterBlockType(toolinit_tools[i].name);
        toolinit_tools[i].prepare(conf);
    }
}

