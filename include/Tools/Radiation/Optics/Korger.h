#ifndef _RADIATION_KORGER_H
#define _RADIATION_KORGER_H

#include "Tools/Radiation/Optics/Korger.h"

namespace __Radiation {
    class Korger {
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
