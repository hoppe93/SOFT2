#ifndef _RADIATION_IMAGE_H
#define _RADIATION_IMAGE_H

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/Output/RadiationOutput.h"
#include "Tools/Radiation/Output/RadiationOutputException.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class Image : public RadiationOutput {
        private:
            std::string output;
            int nrowpixels, ncolpixels, ntotpixels;

            slibreal_t *image, *imageQ, *imageU, *imageV;
            static slibreal_t *global_image,
                *global_imageQ, *global_imageU, *global_imageV;

            bool storeStokesParameters = false;

			static const std::string DEFAULT_QUANTITIES[];
			static const unsigned int NDEFAULT_QUANTITIES;
        public:
            Image(Detector*, MagneticField2D*, ParticleGenerator*);
            ~Image();

            virtual void Configure(ConfigBlock*, ConfigBlock*) override;
            virtual void Finish() override;
            virtual void Handle(Detector*, Model*, RadiationParticle*) override;
            virtual void Initialize() override;
            virtual void Welcome(const std::string&) override;

            // Only called on root thread
            void Generate() override;

            void AllocateImage();
            void GetImagePixel(Detector*, RadiationParticle*, int&, int&);
            static void GetImagePixel(Detector*, RadiationParticle*, int, int, int&, int&);

            virtual bool MeasuresPolarization() override { return storeStokesParameters; }
    };

    class ImageException : public RadiationOutputException {
        public:
            template<typename ... Args>
            ImageException(const std::string &msg, Args&& ... args)
                : RadiationOutputException(msg, std::forward<Args>(args) ...) {
                AddModule("Image");
            }
    };
}

#endif/*_RADIATION_IMAGE_H*/
