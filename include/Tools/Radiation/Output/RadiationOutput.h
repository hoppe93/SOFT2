#ifndef _RADIATION_OUTPUT_H
#define _RADIATION_OUTPUT_H

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
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class RadiationOutput {
        private:
            std::string name;
			std::vector<std::string> common_quantities;

			typedef std::function<void(SFile*)> qfunc;
			std::map<std::string, qfunc> all_quantities;

			// Private functions
			void write_domain(SFile*, const std::string&);
        protected:
            Detector *detector;
            MagneticField2D *magfield;
			ParticleGenerator *particlegenerator;

			static const char
				DETECTOR_APERTURE[],
				DETECTOR_DIRECTION[],
				DETECTOR_EHAT1[],
				DETECTOR_EHAT2[],
				DETECTOR_POSITION[],
				DETECTOR_VISANG[],
				PARAM1[],
				PARAM1NAME[],
				PARAM2[],
				PARAM2NAME[],
				RO_DOMAIN[],
				RO_WALL[],
				R[],
				TP_BOUNDARY[];
        public:
            RadiationOutput(Detector *d, MagneticField2D *m, ParticleGenerator *pgen) {
                this->detector = d;
                this->magfield = m;
				this->particlegenerator = pgen;
				this->InitializeCommonQuantities();
            }
            virtual ~RadiationOutput() { }
            virtual void Configure(ConfigBlock*, ConfigBlock*) = 0;
            virtual void Finish() = 0;
            virtual void Generate() = 0;
            virtual void Handle(Detector*, Model*, RadiationParticle*) = 0;
            virtual void Initialize() = 0;
            virtual void Welcome(const std::string&) = 0;

            virtual bool MeasuresPolarization() = 0;

			void ConfigureCommonQuantities(const std::string*, const unsigned int, const std::vector<std::string>& list=std::vector<std::string>());
			void InitializeCommonQuantities();
			void WriteCommonQuantities(SFile*);

            const std::string &GetName() const { return this->name; }
            void SetName(const std::string &name) { this->name = name; }
    };
}

#endif/*_RADIATION_OUTPUT_H*/
