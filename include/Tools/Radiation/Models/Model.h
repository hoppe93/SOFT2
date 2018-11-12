#ifndef _RADIATION_MODEL_H
#define _RADIATION_MODEL_H

namespace __Radiation {
    class Model;
};

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Orbit/Orbit.h"
#include "SOFT.h"
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class Model {
        private:
        protected:
            Radiation *parent;
            bool nonzero=false;

        public:
            Model(Radiation *rad) { this->parent=rad; }
            ~Model() {}

            virtual void Configure(struct global_settings*, ConfigBlock*, ConfigBlock*) = 0;
            virtual void InitializeOrbit(Orbit*) = 0;
            virtual void InitializeTimestep(RadiationParticle*) = 0;

            virtual void HandleParticle(RadiationParticle*, const slibreal_t, const slibreal_t) = 0;

            virtual const std::string GetDescription() const = 0;

            Radiation *GetParent() { return this->parent; }

            virtual slibreal_t GetPower() = 0;
            virtual slibreal_t GetPowerQ() = 0;
            virtual slibreal_t GetPowerU() = 0;
            virtual slibreal_t GetPowerV() = 0;
            virtual slibreal_t *GetSpectrum() = 0;
            virtual slibreal_t *GetStokesI() = 0;
            virtual slibreal_t *GetStokesQ() = 0;
            virtual slibreal_t *GetStokesU() = 0;
            virtual slibreal_t *GetStokesV() = 0;
            virtual slibreal_t *GetWavelengths() = 0;
            virtual unsigned int GetNWavelengths() = 0;

            bool IsNonZero() { return this->nonzero; }

            /****************
             * MODEL COLORS *
             ****************/
            const std::string CONE_COLOR      = "\e[31;1m"; // Red + bold
            const std::string ANGDIST_COLOR   = "\e[34;1m"; // Blue + bold
            const std::string ISOTROPIC_COLOR = "\e[33;1m"; // Yellow + bold

            /*******************
             * EMISSION COLORS *
             *******************/
            const std::string BREMSSTRAHLUNG_COLOR = "\e[91m"; // Light red
            const std::string SYNCHROTRON_COLOR    = "\e[94m"; // Ligh blue
            const std::string UNIT_EMISSION_COLOR  = "\e[93m"; // Light yellow
    };
};

#endif/*_RADIATION_MODEL_H*/
