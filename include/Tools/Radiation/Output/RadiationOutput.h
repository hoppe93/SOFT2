#ifndef _RADIATION_OUTPUT_H
#define _RADIATION_OUTPUT_H

#include "Tools/OutputModule.h"

namespace __Radiation {
    class RadiationOutput;
}

#include <functional>
#include <map>
#include <string>
#include <vector>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/SFile.h>
#include "Tools/Radiation/Models/Model.h"
#include "PhaseSpace/ParticleGenerator.h"
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class RadiationOutput : public OutputModule {
        protected:
            Detector *detector;

            static const char
                DETECTOR_APERTURE[],
                DETECTOR_DIRECTION[],
                DETECTOR_EHAT1[],
                DETECTOR_EHAT2[],
                DETECTOR_POSITION[],
                DETECTOR_ROLL[],
                DETECTOR_VISANG[];

        public:
            RadiationOutput(
				Detector *d, MagneticField2D *m,
				ParticleGenerator *pgen, SOFT *soft
			) : OutputModule(m, pgen, soft) {
                this->detector = d;

                this->InitializeCommonQuantities();
            }

            virtual void Configure(ConfigBlock*, ConfigBlock*) = 0;
            virtual void Finish() = 0;
            virtual void Generate() = 0;
            virtual void Handle(Detector*, Model*, RadiationParticle*) = 0;
            virtual void Initialize() = 0;
            virtual void Welcome(const std::string&) = 0;

            virtual bool MeasuresPolarization() = 0;

            virtual void InitializeCommonQuantities() override;
    };
}

#endif/*_RADIATION_OUTPUT_H*/
