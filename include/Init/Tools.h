#ifndef _INIT_TOOLS_H
#define _INIT_TOOLS_H

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>

#include "Orbit/ParticlePusher.h"
#include "PhaseSpace/ParticleGenerator.h"
#include "Tools/Tool.h"

struct soft_toolspec {
    std::string name;
    confblock_t id;
    Tool *(*init)(SOFT*, struct global_settings*, ConfigBlock*, ConfigBlock*, ParticleGenerator*, ParticlePusher*, MagneticField2D*);
    void (*prepare)(Configuration*);
};

void InitToolConfig(Configuration*);

#endif/*_INIT_TOOLS_H*/
