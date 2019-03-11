#ifndef _CONE_EMISSION_H
#define _CONE_EMISSION_H

#include <softlib/config.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class ConeEmission {
        protected:
            slibreal_t power, totQ, totU, totV;
            slibreal_t *I, *Q, *U, *V, *wavelengths;
            unsigned int nwavelengths = 0;

            Detector *detector;
            MagneticField2D *magfield;
        public:
            ConeEmission(Detector*, MagneticField2D*);
            virtual ~ConeEmission();
            virtual void HandleParticle(RadiationParticle*, bool) = 0;

            slibreal_t *GetSpectrum() { return this->I; }
            slibreal_t GetTotalEmission() { return this->power; }
            slibreal_t *GetWavelengths() { return this->wavelengths; }
            unsigned int GetNWavelengths() { return this->nwavelengths; }

            slibreal_t GetPowerQ()   { return this->totQ; }
            slibreal_t GetPowerU()   { return this->totU; }
            slibreal_t GetPowerV()   { return this->totV; }
            slibreal_t *GetStokesI() { return this->I; }
            slibreal_t *GetStokesQ() { return this->Q; }
            slibreal_t *GetStokesU() { return this->U; }
            slibreal_t *GetStokesV() { return this->V; }
    };
}

#endif/*_CONE_EMISSION_H*/
