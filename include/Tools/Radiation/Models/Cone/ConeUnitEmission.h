#ifndef _CONE_UNIT_EMISSION_H
#define _CONE_UNIT_EMISSION_H

#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Cone/ConeEmission.h"
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class ConeUnitEmission : public ConeEmission {
        public:
            ConeUnitEmission(Detector *det, MagneticField2D *mf) : ConeEmission(det, mf) {};
            ~ConeUnitEmission() {}

            void HandleParticle(RadiationParticle*, bool);

            void CalculateSpectrum();
            void CalculateTotalEmission();
    };

    class ConeUnitException : public ConeException {
        public:
            template<typename ... Args>
            ConeUnitException(const std::string &msg, Args&& ... args)
                : ConeException(msg, std::forward<Args>(args) ...) {
                AddModule("UnitEmission");
            }
    };
}

#endif/*_CONE_UNIT_EMISSION_H*/
