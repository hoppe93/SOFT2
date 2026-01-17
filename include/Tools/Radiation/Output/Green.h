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
        public:
            enum MPI_Output_Mode {
                CHUNKED,
                CONTIGUOUS
            };
        private:
            std::string output;
            std::string format;

			static const std::string DEFAULT_QUANTITIES[];
			static const unsigned int NDEFAULT_QUANTITIES;

            slibreal_t *global_function;
            slibreal_t *function;
            size_t fsize, fsizeWithoutStokes, *factors;
            sfilesize_t *dimensions, ndimensions;

            slibreal_t *rgrid, *p1grid, *p2grid;
            int p1type, p2type;

            int nrowpixels, ncolpixels,
                suboffseti, suboffsetj,
                subnrowpixels, subncolpixels,
                nformat;
            bool storeStokesParameters;
            bool storeFAsLinearArray=false;
            bool
                hasP1=false, hasP2=false, hasR=false,
                hasI=false, hasJ=false, hasW=false,
                weighWithDistribution = false,
                withJacobian = true,
                containsAllPhaseSpaceParameters = false,
                pixelsset = false;      // True of the number of pixels has been set

            // Indices of parameters
            int
                i1, i2, ii, ij, ir, iw;

            int
                nr, n1, n2, nw;
            
            unsigned int
                end_r, end_1, end_2,
                start_r, start_1, start_2;
            
            enum ParticleGenerator::MPI_Distribute_Mode mpi_distribute_mode;
            enum MPI_Output_Mode mpi_output_mode = MPI_Output_Mode::CONTIGUOUS;
        public:
            Green(Detector*, MagneticField2D*, ParticleGenerator*, SOFT*);
            ~Green();

            virtual void Configure(ConfigBlock*, ConfigBlock*) override;
            virtual void Finish() override;
            virtual void Handle(Detector*, Model*, RadiationParticle*) override;
            virtual void Initialize() override;
            virtual void Welcome(const std::string&) override;

            // Only called on root thread
            void Generate() override;

            std::string GetChunkedName(const int);

            void PrepareAllocateGreen();
            bool GetIndex(__Radiation::Detector*, __Radiation::RadiationParticle*, size_t*, size_t*);
            bool MeasuresPolarization() override { return storeStokesParameters; }

            std::string TranslateFormat(const std::string&);
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
