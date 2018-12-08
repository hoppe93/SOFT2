#ifndef _RADIATION_SPECTRUM_H
#define _RADIATION_SPECTRUM_H

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/Output/RadiationOutput.h"
#include "Tools/Radiation/Output/RadiationOutputException.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class Spectrum : public RadiationOutput {
        private:
            std::string output;
            unsigned int nwavelengths;

            slibreal_t *I, *Q, *U, *V, *wavelengths;
            static slibreal_t *global_I, *global_Q, *global_U, *global_V;
        public:
            Spectrum(Detector*, MagneticField2D*, ParticleGenerator*);
            ~Spectrum();

            void Configure(ConfigBlock*, ConfigBlock*);
            void Finish();
            void Handle(Detector*, Model*, RadiationParticle*);

            // Only called on root thread
            void Generate();

            bool MeasuresPolarization() { return false; }
    };

    class SpectrumException : public RadiationOutputException {
        public:
            template<typename ... Args>
            SpectrumException(const std::string &msg, Args&& ... args)
                : RadiationOutputException(msg, std::forward<Args>(args) ...) {
                AddModule("Spectrum");
            }
    };
}

#endif/*_RADIATION_SPECTRUM_H*/
