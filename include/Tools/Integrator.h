#ifndef _SOFT_INTEGRATOR_H
#define _SOFT_INTEGRATOR_H

#include <string>
#include "Tools/OutputModule.h"
#include "Tools/Tool.h"

namespace __SOFT {
    class Integrator : public Tool, public OutputModule {
    private:
        slibreal_t I;
		SOFT *soft;

		std::string output;
    public:
        Integrator(SOFT*, MagneticField2D*, ParticleGenerator*, ParticlePusher*);

        virtual void Configure(struct global_settings*, ConfigBlock*, ConfigBlock*) override;
        virtual void Finish() override {}
        virtual void Handle(Orbit*, Particle*) override;
        virtual void Initialize() override {}
        virtual void Output() override;
        virtual void Welcome(const std::string& prefix="  ") override;

        static void PrepareConfiguration(Configuration*);
    };
}

#endif/*_INTEGRATOR_H*/
