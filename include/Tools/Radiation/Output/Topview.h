#ifndef _RADIATION_TOPVIEW_H
#define _RADIATION_TOPVIEW_H

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/Output/RadiationOutput.h"
#include "Tools/Radiation/Output/RadiationOutputException.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class Topview : public RadiationOutput {
        private:
            std::string output;
            unsigned int npixels, ntotpixels;
            slibreal_t max_radius;

            slibreal_t *topview;
            static slibreal_t *global_topview;
        public:
            Topview(Detector*, MagneticField2D*, ParticleGenerator*);
            ~Topview();

            void Configure(ConfigBlock*, ConfigBlock*);
            void Finish();
            void Handle(Detector*, Model*, RadiationParticle*);

            // Only called on root thread
            void Generate();

            void AllocateTopview();
            void GetTopviewPixel(RadiationParticle*, unsigned int*, unsigned int*);

            bool MeasuresPolarization() { return false; }
    };

    class TopviewException : public RadiationOutputException {
        public:
            template<typename ... Args>
            TopviewException(const std::string &msg, Args&& ... args)
                : RadiationOutputException(msg, std::forward<Args>(args) ...) {
                AddModule("Topview");
            }
    };
}

#endif/*_RADIATION_TOPVIEW_H*/
