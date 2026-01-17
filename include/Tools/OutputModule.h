#ifndef _OUTPUT_MODULE_H
#define _OUTPUT_MODULE_H

#include <functional>
#include <map>
#include <string>
#include <vector>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/DistributionFunction/DistributionFunction.h>
#include <softlib/SFile.h>
//#include "Tools/Radiation/Models/Model.h"
#include "PhaseSpace/ParticleGenerator.h"
//#include "Tools/Radiation/RadiationParticle.h"

class OutputModule {
    private:
        std::string name;
        std::vector<std::string> common_quantities;

        typedef std::function<void(SFile*)> qfunc;
        std::map<std::string, qfunc> all_quantities;

        // Private functions
        void write_domain(SFile*, const std::string&);
        void write_separatrix(SFile*, const std::string&);
    protected:
        MagneticField2D *magfield;
        ParticleGenerator *particlegenerator;
		SOFT *soft;

        static const char
            DATETIME[],
			DISTFUNC[],
            PARAM1[],
            PARAM1NAME[],
            PARAM2[],
            PARAM2NAME[],
            RO_DOMAIN[],
            RO_SEPARATRIX[],
            RO_WALL[],
            R[],
            VERSION_SOFT[];
            //VERSION_SOFTLIB[];
    public:
        OutputModule(
			MagneticField2D *m, ParticleGenerator *pgen,
			SOFT *soft
		) {
            this->magfield = m;
            this->particlegenerator = pgen;
			this->soft = soft;
            this->InitializeCommonQuantities();
        }
        virtual ~OutputModule() { }

        virtual void ConfigureCommonQuantities(const std::string*, const unsigned int, const std::vector<std::string>& list=std::vector<std::string>());
        void DefineCommonQuantity(const std::string&, const qfunc&);
        virtual void InitializeCommonQuantities();
        void WriteCommonQuantities(SFile*);

        const std::string &GetName() const { return this->name; }
        void SetName(const std::string &name) { this->name = name; }
};

#endif/*_OUTPUT_MODULE_H*/
