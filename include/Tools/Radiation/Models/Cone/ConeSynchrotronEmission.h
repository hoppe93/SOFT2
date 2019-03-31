#ifndef _CONE_SYNCHROTRON_EMISSION_H
#define _CONE_SYNCHROTRON_EMISSION_H

#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Cone/ConeEmission.h"
#include "Tools/Radiation/Optics/Optics.h"
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class ConeSynchrotronEmission : public ConeEmission {
        private:
            struct Optics::Efield Efield;
        public:
            ConeSynchrotronEmission(Detector*, MagneticField2D*, bool polarization=false);
            ~ConeSynchrotronEmission();

            void HandleParticle(RadiationParticle*, bool);

            void CalculatePolarization(RadiationParticle*);
            void CalculateSpectrum(RadiationParticle*);
            void CalculateTotalEmission(RadiationParticle*);

            void IntegrateSpectrum();
            void IntegrateSpectrumStokes();
    };

    class ConeSynchrotronException : public ConeException {
        public:
            template<typename ... Args>
            ConeSynchrotronException(const std::string &msg, Args&& ... args)
                : ConeException(msg, std::forward<Args>(args) ...) {
                AddModule("Synchrotron");
            }
    };
}

#endif/*_CONE_SYNCHROTRON_EMISSION_H*/
