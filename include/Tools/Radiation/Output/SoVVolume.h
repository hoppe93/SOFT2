#ifndef _RADIATION_SOVVOLUME_H
#define _RADIATION_SOVVOLUME_H

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/Output/RadiationOutput.h"
#include "Tools/Radiation/Output/RadiationOutputException.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class SoVVolume : public RadiationOutput {
        private:
            std::string output;
            size_t np1, np2;
            slibreal_t *p1grid, *p2grid;

            int p1type, p2type;

            slibreal_t *volumearray=nullptr;
            size_t arraysize, narrayelements;

            static slibreal_t *global_volumearray;

			static const std::string DEFAULT_QUANTITIES[];
			static const unsigned int NDEFAULT_QUANTITIES;

        public:
            SoVVolume(Detector*, MagneticField2D*, ParticleGenerator*, SOFT*);
            virtual ~SoVVolume();

            virtual void Configure(ConfigBlock*, ConfigBlock*) override;
            virtual void Finish() override;
            virtual void Generate() override;
            virtual void Handle(Detector*, Model*, RadiationParticle*) override;
            virtual void Initialize() override;
            virtual void Welcome(const std::string&) override;

            void AllocateVolume();

            virtual bool MeasuresPolarization() override { return false; }
    };

    class SoVVolumeException : public RadiationOutputException {
        public:
            template<typename ... Args>
            SoVVolumeException(const std::string &msg, Args&& ... args)
                : RadiationOutputException(msg, std::forward<Args>(args) ...) {
                AddModule("SoVVolume");
            }
    };
}

#endif/*_RADIATION_SOVVOLUME_H*/
