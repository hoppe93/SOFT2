#ifndef _CONE_SYNCHROTRON_EMISSION_H
#define _CONE_SYNCHROTRON_EMISSION_H

#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Cone/ConeEmission.h"
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

slibreal_t synchrotron_func1(const slibreal_t);
slibreal_t synchrotron_func2(const slibreal_t);

namespace __Radiation {
    class ConeSynchrotronEmission : public ConeEmission {
        protected:
            template<bool calculatePolarization>
            void __CalculateSpectrum(RadiationParticle*);

        public:
            ConeSynchrotronEmission(Detector *det) : ConeEmission(det) {};
            ~ConeSynchrotronEmission() {}

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
