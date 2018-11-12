#ifndef _CONE_EMISSION_H
#define _CONE_EMISSION_H

#include <softlib/config.h>
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class ConeEmission {
        protected:
            slibreal_t power;
            slibreal_t *I, *Q, *U, *V, *wavelengths;
            unsigned int nwavelengths = 0;

            Detector *detector;
        public:
            ConeEmission(Detector*);
            virtual ~ConeEmission();
            virtual void HandleParticle(RadiationParticle*, bool) = 0;

            slibreal_t *GetSpectrum() { return this->I; }
            slibreal_t GetTotalEmission() { return this->power; }
            slibreal_t *GetWavelengths() { return this->wavelengths; }
            unsigned int GetNWavelengths() { return this->nwavelengths; }

            slibreal_t *GetStokesI() { return this->I; }
            slibreal_t *GetStokesQ() { return this->Q; }
            slibreal_t *GetStokesU() { return this->U; }
            slibreal_t *GetStokesV() { return this->V; }
    };
}

#endif/*_CONE_EMISSION_H*/
