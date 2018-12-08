#ifndef _RADIATION_GREEN_H
#define _RADIATION_GREEN_H

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "PhaseSpace/ParticleGenerator.h"
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/Output/RadiationOutput.h"
#include "Tools/Radiation/Output/RadiationOutputException.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class Green : public RadiationOutput {
        private:
            std::string output;
            std::string format;

            static slibreal_t *global_function;
            slibreal_t *function;
            size_t fsize, fsizeWithoutStokes, *factors;

            slibreal_t *rgrid, *p1grid, *p2grid;
            int p1type, p2type;

            int nrowpixels, ncolpixels,
                suboffseti, suboffsetj,
                subnrowpixels, subncolpixels,
                nformat;
            bool storeStokesParameters;
            bool
                hasP1=false, hasP2=false, hasR=false,
                hasI=false, hasJ=false, hasW=false,
                weighWithDistribution = false,
                withJacobian = true,
                containsAllMomentumSpaceParameters = false,
                pixelsset = false;      // True of the number of pixels has been set

            // Indices of parameters
            int
                i1, i2, ii, ij, ir, iw;

            int
                nr, n1, n2, nw;
        public:
            Green(Detector*, MagneticField2D*, ParticleGenerator*);
            ~Green();

            virtual void Configure(ConfigBlock*, ConfigBlock*) override;
            virtual void Finish() override;
            virtual void Handle(Detector*, Model*, RadiationParticle*) override;
            virtual void Initialize() override;
            virtual void Welcome(const std::string&) override;

            // Only called on root thread
            void Generate();

            void PrepareAllocateGreen();
            void GetIndex(__Radiation::Detector*, __Radiation::RadiationParticle*, size_t*, size_t*);
            bool MeasuresPolarization() { return storeStokesParameters; }

            void ValidateFormat(const std::string&);
            void ValidateSubPixels();
    };

    class GreenException : public RadiationOutputException {
        public:
            template<typename ... Args>
            GreenException(const std::string &msg, Args&& ... args)
                : RadiationOutputException(msg, std::forward<Args>(args) ...) {
                AddModule("Green");
            }
    };
}

#endif/*_RADIATION_GREEN_H*/
