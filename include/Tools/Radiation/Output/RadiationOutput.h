#ifndef _RADIATION_OUTPUT_H
#define _RADIATION_OUTPUT_H

namespace __Radiation {
    class RadiationOutput;
}

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Tools/Radiation/Models/Model.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class RadiationOutput {
        private:
            std::string name;
        protected:
            Detector *detector;
            MagneticField2D *magfield;
        public:
            RadiationOutput(Detector *d, MagneticField2D *m) {
                this->detector = d;
                this->magfield = m;
            }
            virtual ~RadiationOutput() { }
            virtual void Configure(ConfigBlock*, ConfigBlock*) = 0;
            virtual void Finish() = 0;
            virtual void Generate() = 0;
            virtual void Handle(Detector*, Model*, RadiationParticle*) = 0;
            virtual void Initialize() = 0;
            virtual void Welcome(const std::string&) = 0;

            virtual bool MeasuresPolarization() = 0;

            const std::string &GetName() const { return this->name; }
            void SetName(const std::string &name) { this->name = name; }
    };
}

#endif/*_RADIATION_OUTPUT_H*/
