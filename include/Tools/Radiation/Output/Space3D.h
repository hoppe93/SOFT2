#ifndef _RADIATION_SPACE3D_H
#define _RADIATION_SPACE3D_H

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/Output/RadiationOutput.h"
#include "Tools/Radiation/Output/RadiationOutputException.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class Space3D : public RadiationOutput {
        private:
            std::string output;
            Vector<3> point0, point1;
            long long int pixelsX, pixelsY, pixelsZ;

            slibreal_t *s3dimage=nullptr;
            size_t imagesize;

			static const std::string DEFAULT_QUANTITIES[];
			static const unsigned int NDEFAULT_QUANTITIES;

        public:
            Space3D(Detector*, MagneticField2D*, ParticleGenerator*);
            virtual ~Space3D();

            virtual void Configure(ConfigBlock*, ConfigBlock*) override;
            virtual void Finish() override;
            virtual void Generate() override;
            virtual void Handle(Detector*, Model*, RadiationParticle*) override;
            virtual void Initialize() override;
            virtual void Welcome(const std::string&) override;

            void GetPoint(ConfigBlock*, const std::string&, Vector<3>&);

            virtual bool MeasuresPolarization() override { return false; }
    };

    class Space3DException : public RadiationOutputException {
        public:
            template<typename ... Args>
            Space3DException(const std::string &msg, Args&& ... args)
                : RadiationOutputException(msg, std::forward<Args>(args) ...) {
                AddModule("Space3D");
            }
    };
}

#endif/*_RADIATION_SPACE3D_H*/
