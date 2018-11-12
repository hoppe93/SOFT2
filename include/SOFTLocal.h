#ifndef _SOFT_LOCAL_H
#define _SOFT_LOCAL_H

#include <vector>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/DistributionFunction/DistributionFunction.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "PhaseSpace/ParticleGenerator.h"
#include "Orbit/ParticlePusher.h"
#include "Tools/ToolHandler.h"
#include "SOFT.h"

/**
 * A SOFT object designed to run asynchronously
 * with other SOFT objects.
 */

class SOFTLocal {
    private:
        unsigned int id;        // Object ID (typically thread number)
    public:
        // Pointers to global objects
        ParticleGenerator *partgen;
        SOFT *soft;
        DistributionFunction *distribution;

        // Local objects
        MagneticField2D *magfield;
        ParticlePusher *pusher;
        ToolHandler *thandler;
        
        SOFTLocal(SOFT*, unsigned int, bool clone=true);

        unsigned int GetID() { return this->id; }

        // Initialization of sub-modules
        static ParticlePusher *InitPusher(MagneticField2D*, struct global_settings*, Configuration*);
        static ToolHandler *InitTools(struct global_settings*, Configuration*, ParticleGenerator*, ParticlePusher*, MagneticField2D*);

        void Finish();
        void Run(unsigned int *invalid=nullptr);
        void Output();
        void Welcome();
};

#endif/*_SOFT_LOCAL_H*/
