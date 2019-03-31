#ifndef _RADIATION_KORGER_H
#define _RADIATION_KORGER_H

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/Optics/Optics.h"

namespace __Radiation {
    class Korger : public Optics {
        public:
            Korger(Detector*, ConfigBlock*);

            virtual void ApplyOptics(
                const struct Optics::Efield&,
                slibreal_t *I, slibreal_t *Q=nullptr,
                slibreal_t *U=nullptr, slibreal_t *V=nullptr
            ) override;

            void Configure(ConfigBlock*);
    };
}

#endif/*_RADIATION_KORGER_H*/
