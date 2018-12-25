#ifndef _ANGULAR_DISTRIBUTION_EMISSION_H
#define _ANGULAR_DISTRIBUTION_EMISSION_H

#include <softlib/config.h>
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationException.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class ADEmission {
        protected:
            slibreal_t power, powerQ, powerU, powerV;
            slibreal_t *I, *Q, *U, *V, *wavelengths;
            unsigned int nwavelengths = 0;
        public:
            ADEmission(Detector*);
            virtual ~ADEmission();
            virtual slibreal_t Evaluate(Vector<3>&, slibreal_t, slibreal_t, bool) = 0;
            virtual void InitializeToroidalStep(const slibreal_t, const slibreal_t) = 0;
            virtual void Prepare(RadiationParticle*, bool) = 0;

            virtual void CalculateAngularDistribution(Vector<3>&, slibreal_t, slibreal_t) = 0;
            virtual void CalculatePolarization(Vector<3>&, slibreal_t, slibreal_t) = 0;
            virtual void CalculateSpectrum(Vector<3>&, slibreal_t, slibreal_t) = 0;

            slibreal_t *GetSpectrum() const { return this->I; }
            slibreal_t GetTotalEmission() const { return this->power; }
            slibreal_t GetTotalEmissionQ() const { return this->powerQ; }
            slibreal_t GetTotalEmissionU() const { return this->powerU; }
            slibreal_t GetTotalEmissionV() const { return this->powerV; }
            slibreal_t *GetWavelengths() const { return this->wavelengths; }
            unsigned int GetNWavelengths() const { return nwavelengths; }

            slibreal_t *GetStokesI() const { return this->I; }
            slibreal_t *GetStokesQ() const { return this->Q; }
            slibreal_t *GetStokesU() const { return this->U; }
            slibreal_t *GetStokesV() const { return this->V; }
    };
}

#endif/*_ANGULAR_DISTRIBUTION_EMISSION_H*/
