/**
 * Initialize the particle pusher object.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>

#include "Init/InitConfig.h"
#include "Orbit/ParticlePusher.h"
#include "SOFTLocal.h"
#include "SOFTException.h"

/**
 * Initialize the particle pusher.
 *
 * magfield: Magnetic field object to use for initialization.
 * globset:  Global SOFT settings.
 * input:    Configuration object to initialize from.
 */
ParticlePusher *SOFTLocal::InitPusher(MagneticField2D *magfield, struct global_settings *globset, Configuration *input) {
    ConfigBlock *root = input->GetRootBlock();
    ConfigBlock *cfb = root->GetConfigBlock(CONFBLOCK_PARTICLEPUSHER, globset->particle_pusher);
    if (cfb == nullptr)
        throw SOFTException("Particle pusher '%s' not defined.", globset->particle_pusher.c_str());

    return new ParticlePusher(magfield, globset, cfb, root);
}

