/**
 * Implementation of a general integrator.
 */

#include "Tools/Integrator.h"


/**
 * Constructor.
 */
__SOFT::Integrator::Integrator(SOFT *soft, MagneticField2D *mf, ParticleGenerator *pgen, ParticlePusher *pp)
    : Tool("Integrator"), OutputModule(mf, pgen, soft), soft(soft) {
    
    this->I = 0;
}

/**
 * Configure this module.
 */
void __SOFT::Integrator::Configure(
	struct global_settings*, ConfigBlock *conf, ConfigBlock*
) {
	this->Tool::SetName(conf->GetName());
	this->OutputModule::SetName(conf->GetName());

	if (conf->HasSetting("output"))
		this->output = (*conf)["output"];
}

/**
 * Handle the given particle orbit.
 */
void __SOFT::Integrator::Handle(Orbit *o, Particle *p) {
    unsigned int ntau = o->GetNTau();
    slibreal_t s = 0;
    const slibreal_t *Jdtdrho = o->GetJdtdrho();

    // Evaluate time integral (<=> poloidal integral)
    for (unsigned int i = 0; i < ntau; i++)
        s += Jdtdrho[i];

    slibreal_t
        dp1   = p->GetJMomentum1(),
        dp2   = p->GetJMomentum1();
    unsigned int
        np1   = p->GetNP1(),
        np2   = p->GetNP2(),
        nzeta = p->GetNZeta();

    // Only include gyro angle integral if we're
    // integrating over momentum
    slibreal_t Ig = (np1==1 && np2==1 && nzeta==1 ? 1 : 2*M_PI);

    s *=
        2*M_PI*             // Toroidal integral
        Ig*                 // Gyro angle integral
        p->GetF()*          // Distribution function
        p->GetJMomentum1()* // Momentum coordinate 1 jacobian & differential element (J1*dp1)
        p->GetJMomentum2(); // Momentum coordinate 2 jacobian & differential element (J2*dp2)

    #pragma omp atomic
    this->I += s;
}

/**
 * Generate output from this module.
 */
void __SOFT::Integrator::Output() {
	if (this->output.empty())
		SOFT::PrintInfo("Integral:  %.16e", this->I);
	else {
		SFile *sf = SFile::Create(this->output, SFILE_MODE_WRITE);

		sf->WriteScalar("V", this->I);

		sf->Close();
	}
}

/**
 * Called when constructing a new 'Integrator' object.
 */
void __SOFT::Integrator::PrepareConfiguration(Configuration*) {}

