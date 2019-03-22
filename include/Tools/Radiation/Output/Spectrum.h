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

			static const std::string DEFAULT_QUANTITIES[];
			static const unsigned int NDEFAULT_QUANTITIES;

        public:
            Spectrum(Detector*, MagneticField2D*, ParticleGenerator*);
            ~Spectrum();

            virtual void Configure(ConfigBlock*, ConfigBlock*) override;
            virtual void Finish() override;
            virtual void Handle(Detector*, Model*, RadiationParticle*) override;
            virtual void Initialize() override;
            virtual void Welcome(const std::string&) override;

            // Only called on root thread
            virtual void Generate() override;

            virtual bool MeasuresPolarization() override { return false; }
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
