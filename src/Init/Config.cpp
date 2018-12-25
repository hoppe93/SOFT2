/**
 * Initialize configuration object.
 */

#include <softlib/Configuration.h>
#include "Init/Init.h"
#include "Init/InitConfig.h"
#include "Init/Tools.h"

int CONFBLOCK_EQUATION_GC,
    CONFBLOCK_EQUATION_PARTICLE,
    CONFBLOCK_MAGNETICFIELD,
    CONFBLOCK_PARTICLEGENERATOR,
    CONFBLOCK_PARTICLEPUSHER,
    CONFBLOCK_DISTRIBUTION,
    CONFBLOCK_RADIALPROFILE;

/**
 * Initializes the given configuration object.
 * This function should be called before any configuration
 * file has been loaded so that the appropriate blocks
 * can be registered.
 *
 * conf: Freshly created configuration object, which
 *       should _not_ have been initialized with a
 *       configuration file yet.
 */
void InitConfig(Configuration *conf) {
    CONFBLOCK_EQUATION_GC       = conf->RegisterBlockType("@EquationGuidingCenter");
    CONFBLOCK_EQUATION_PARTICLE = conf->RegisterBlockType("@EquationParticle");
    CONFBLOCK_MAGNETICFIELD     = conf->RegisterBlockType("@MagneticField");
    CONFBLOCK_PARTICLEGENERATOR = conf->RegisterBlockType("@ParticleGenerator");
    CONFBLOCK_PARTICLEPUSHER    = conf->RegisterBlockType("@ParticlePusher");
    CONFBLOCK_DISTRIBUTION      = conf->RegisterBlockType("@DistributionFunction");
    CONFBLOCK_RADIALPROFILE     = conf->RegisterBlockType("@RadialProfile");

    InitToolConfig(conf);
}

