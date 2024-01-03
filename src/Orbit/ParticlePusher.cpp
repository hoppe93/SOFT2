/**
 * ParticlePusher
 * High-level interface for pushing particles/guiding-centers.
 *
 * This class encapsulates all logic for generating particle/guiding-center
 * orbits.
 */

#include <string>

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/Configuration/ConfigurationScript.h>
#include <softlib/Integrator.h>
#include <softlib/IntegratorEquation.h>
#include <softlib/RKDP45.h>

#include "Init/Init.h"
#include "Init/InitConfig.h"
#include "Orbit/GuidingCenterEquation.h"
#include "Orbit/ParticleEquation.h"
#include "Orbit/ParticlePusher.h"
#include "PhaseSpace/Particle.h"

using namespace std;

/**
 * DEFAULT CONFIGURATION
 * This string defines the default configuration
 * for this module. Note that also unused blocks
 * can (or rather should) be defined here to
 * provide default settings for all possible
 * options.
 */
const string ParticlePusher_config =
"equation=guiding-center;\n"
"timeunit=poloidal;\n"
"time=1;\n"
"nt=1e3;\n"
"force_numerical_jacobian=no;\n"
"nudgevalue=__default__;\n";

const string ParticlePusher::equation_defaults =
"@Equation guiding-center (guiding-center) {\n"
"   method=rkdp45;\n"
"   tolerance=1e-9;\n"
"}\n"
"@Equation particle (particle) {\n"
"   method=rkdp45;\n"
"   tolerance=1e-7;\n"
"}\n"
;

/**
 * Constructor.
 *
 * mf:      Magnetic field with which to initialize the new ParticlePusher.
 * globset: Global SOFT settings.
 * conf:    Configuration block containing settings for the ParticlePusher.
 * eqnconf: Configuration (root) block containing all '@Equation' blocks
 *          that equations can be initialized from.
 */
ParticlePusher::ParticlePusher(MagneticField2D *mf, struct global_settings *globset)
    : ParticlePusher(mf, globset, (ConfigBlock*)nullptr, (ConfigBlock*)nullptr) {}
ParticlePusher::ParticlePusher(
    MagneticField2D *mf, struct global_settings *globset,
    ConfigBlock *conf, ConfigBlock *eqnconf
) {
    this->magfield = mf;
    this->globset = globset;
    ResetDomainCheck();
    ResetPoloidalTime();

	InitDefaults();

	// Handle settings
	if (conf != nullptr) {
		vector<string> *v = settings->Merge(conf);
		if (v != nullptr) {
			int i, l = v->size();
			for (i = 0; i < l; i++)
				SOFT::PrintWarning(SOFT::WARNING_OPP_UNRECOGNIZED_SETTING, "ParticlePusher: Unrecognized setting '%s' ignored.", (*v)[i].c_str());
		}
	}

	// Interpret settings
	InitEquation((*settings)["equation"], *eqnconf);

    // Override nude value
    if ((*settings)["nudgevalue"] == "__default__") {
        // Set the 'nudge' value (radial distance between
        // orbits when calculating Jacobians).
        // Note that 'integrator_tol' is initialized when
        // the integrator is set, in 'InitEquation'.
        this->nudge_value = sqrt(this->integrator_tol * mf->GetMagneticAxisR());
    } else {
        Setting *s = settings->GetSetting("nudgevalue");
        if (!s->IsScalar())
            throw ParticlePusherException("Invalid value assigned to parameter 'nudgevalue'. Expected real scalar.");

        this->nudge_value = s->GetScalar();
    }

    // Force numerical calculation of Jacobian determinant?
    Setting *s = settings->GetSetting("force_numerical_jacobian");
    if (!s->IsBool())
        throw ParticlePusherException("Invalid value assigned to parameter 'force_numerical_jacobian'. Expected boolean value.");
    else
        this->forceNumericalJacobian = s->GetBool();

    // Max time
    this->maxtime = init_get_scalar(settings, "time", "ParticlePusher");

    // Time unit
    if ((*settings)["timeunit"] == "poloidal") {
        if (this->equation->GetOrbitType() != ORBIT_TYPE_GUIDING_CENTER) {
            // Set up "timing equation"
            SetupTimingIntegrator();
            this->timeunit = ORBITTIMEUNIT_SECONDS;
        } else
            this->timeunit = ORBITTIMEUNIT_POLOIDAL;
    } else if ((*settings)["timeunit"] == "seconds")
        this->timeunit = ORBITTIMEUNIT_SECONDS;
    else
        throw ParticlePusherException("Unrecognized time unit: '%s'.", (*settings)["timeunit"].c_str());

    // Set domain bounds
    mf->GetDomainBounds(domain_rmin, domain_rmax, domain_zmin, domain_zmax);

    // Number of time steps
    this->ntimesteps = init_get_uint32(settings, "nt", "ParticlePusher");
    this->retorbit = new Orbit(this->equation->GetOrbitType(), this->ntimesteps, globset->include_drifts);
    this->solution = new slibreal_t[this->ntimesteps];
}

/**********************************
 *        INITIALIZATION          *
 **********************************/
/**
 * Create a settings object and populate
 * it with default values.
 */
void ParticlePusher::InitDefaults() {
    ConfigurationScript c;
	settings = c.FromString(ParticlePusher_config);
}

/**
 * Initialize the method to use for
 * particle pushing.
 *
 * equation: Name of equation to use.
 * eqnconf:  ConfigBlock containing '@Equation' blocks
 *           to be used for configuration.
 */
void ParticlePusher::InitEquation(const string& equation, ConfigBlock& eqnconf) {
	ConfigBlock *conf;

    if (eqnconf.HasSubBlock(CONFBLOCK_EQUATION, equation))
		conf = eqnconf.GetConfigBlock(CONFBLOCK_EQUATION, equation);
    else
		throw ParticlePusherException("No equation named '"+equation+"' specified in configuration.");

    if (conf->GetSecondaryType() == "guiding-center") {
		GuidingCenterEquation *eq = new GuidingCenterEquation(this->magfield, this->globset);
        this->equation = eq;

		// Choose integrator
		this->InitGeneralIntegrator(*conf, eq);
    } else if (conf->GetSecondaryType() == "particle") {
		conf = eqnconf.GetConfigBlock(CONFBLOCK_EQUATION, equation);
		
		ParticleEquation *eq = new ParticleEquation(this->magfield, this->globset);
        this->equation = eq;

		// Choose integrator
		this->InitGeneralIntegrator(*conf, eq);
    } else
        throw ParticlePusherException("Unrecognized equation type '%s'.", conf->GetSecondaryType().c_str());

    // Emit a warning if the particle equations of motion are solved,
    // while 'include_drifts = false'. Since drifts cannot be turned off
    // in the particle equations of motion, SOFT must be run with
    // 'include_drifts = true'.
    if (this->equation->GetOrbitType() == ORBIT_TYPE_PARTICLE && !this->globset->include_drifts)
        SOFT::PrintWarning(
            SOFT::WARNING_OPP_PARTICLE_EQUATION_WO_DRIFTS,
            "Solving the particle equations of motion, but 'include_drifts = false'. This is will lead an inconsistent calculation."
        );
}

/**
 * Set up the 'timing' equation to use along with the full-orbit solver.
 * (This equation is used to determine the "poloidal time" of a particle).
 * 
 */
void ParticlePusher::SetupTimingIntegrator() {
    // Set equation
    GuidingCenterEquation *eq = new GuidingCenterEquation(this->magfield, this->globset);
    this->timingEquation = eq;

    // Set up integrator
    slibreal_t integrator_tol = 1e-6;
    RKDP45<6> *rkdp45 = new RKDP45<6>(integrator_tol);
    rkdp45->SetEquation(eq);

    this->timingIntegrator = rkdp45;
    this->useTimingIntegrator = true;
}

/**
 * Initialize one of the generalized integrators
 * for the particle pusher. Currently, the following
 * algorithms are available:
 *
 *   RKDP45  -- Runge-Kutta of order 4(5) with Dormand-Prince coefficients.
 *
 * conf: Settings for the integrator to initialize.
 * eq:   Pointer to the integrator equation to be integrated.
 */
void ParticlePusher::InitGeneralIntegrator(ConfigBlock& conf, IntegratorEquation<6> *eq) {
    Setting *s;
	if (conf["method"] == "rkdp45") {
		s = conf.GetSetting("tolerance");
		if (!s->IsScalar())
			throw ParticlePusherException("Assigned tolerance value is not a scalar.");

        integrator_tol = s->GetScalar();
        
        if (integrator_tol >= 1.0)
            throw ParticlePusherException("Assigned relative tolerance is greater than one.");
        else if (integrator_tol <= REAL_EPSILON)
            throw ParticlePusherException("Assigned relative tolerance is smaller than the machine epsilon (%e).", REAL_EPSILON);

		RKDP45<6> *rkdp45_1 = new RKDP45<6>(integrator_tol),
                  *rkdp45_2 = new RKDP45<6>(integrator_tol);
		rkdp45_1->SetEquation(eq);
        rkdp45_2->SetEquation(eq);

		this->integrator1 = rkdp45_1;
        this->integrator2 = rkdp45_2;
	} else
		throw ParticlePusherException("Unrecognized integration method '"+conf["method"]+"'.");
}

/**
 * Given an finished integrator object, use the secant method
 * to find the point Z = 0 between the last two timesteps.
 *
 * eqn:  Equation that was solved.
 * intg: Integrator object containing the solved orbit.
 */
slibreal_t ParticlePusher::FindPoloidalTime(SOFTEquation *eqn, Integrator<6> *intg) {
    const unsigned int NVARIABLES = 6;
    unsigned int nt = intg->StepsTaken();
    slibreal_t
        a = intg->TimeAt(nt-2),
        b = intg->TimeAt(nt-1),
        Rm = magfield->GetMagneticAxisR(),
        x, za, zx, f[2*6], t[2*6];

    auto gz = [this,&eqn,&intg,&f,&t](slibreal_t z) {
        intg->OutputDense(1, z, z, f, t);
        return (this->zinit - eqn->GetPositionZ(f));
    };

    // Do one iteration manually...
    slibreal_t *sol = intg->SolutionAt(nt-2);
    for (unsigned int i = 0; i < 2*NVARIABLES; i++)
        f[i] = sol[i];

    za = zinit - eqn->GetPositionZ(f);
    zx = zinit - eqn->GetPositionZ(f+NVARIABLES);

    if (zx*za > 0)
        throw ParticlePusherException(
            "Unable to determine poloidal time for orbit that does not close on itself."
        );

    do {
        x = 0.5*(a+b);
        zx = gz(x);

        if (zx == 0)
            return x;
        else if (zx*za > 0)
            a = x;
        else
            b = x;
    } while (fabs(zx/Rm) > integrator_tol);

    return x;
}

/**
 * Calculates the current time in units of the
 * poloidal orbit time. Note that this method
 * only returns integers, and rounds to the closest
 * nearest smaller value.
 *
 * This method counts the number of times the
 * particle passes z = z0 and returns half that
 * value.
 *
 * r: Particle radial position.
 * z: Particle vertical position.
 */
slibreal_t ParticlePusher::GetPoloidalTime(slibreal_t z) {
    // If the particle doesn't move vertically,
    // poloidal time is not defined
    if (fabs(z-zprev) < integrator_tol) {
        restflag++;
        if (restflag >= 2) {
            return (++zpass);
        }
    } else restflag = 0;

    if ((zprev < zinit && zinit < z) ||
        (zprev > zinit && zinit > z))
        zpass++;

    zprev = z;

    if (zpass == 0)
        return 0.0;
    else
        return (zpass/2);
}

/**
 * Resets the poloidal time counter.
 */
void ParticlePusher::ResetPoloidalTime() {
    this->zpass = 0;
    this->restflag = 0;
}

/**
 * Reset the domain check.
 */
void ParticlePusher::ResetDomainCheck() {
    this->outside_domain_flag = false;
}

/**
 * Returns the radial position of the particle
 * in the current state.
 */
slibreal_t ParticlePusher::GetPositionR(SOFTEquation *eqn, Integrator<6> *intg) {
    return eqn->GetPositionR(
        intg->LastSolution()
    );
}
/**
 * Returns the vertical position of the particle
 * in the current state.
 */
slibreal_t ParticlePusher::GetPositionZ(SOFTEquation *eqn, Integrator<6> *intg) {
    return eqn->GetPositionZ(
        intg->LastSolution()
    );
}

/**
 * Check if the given point lies inside the 'domain bounds'.
 * This is a weaker condition than checking whether the point
 * lies within the 'domain'. The 'domain bound' is basically
 * a bounding box, which contains the domain.
 *
 * R: Radial coordinate of point.
 * Z: Vertical coordinate of point.
 */
bool ParticlePusher::InsideDomainBounds(const slibreal_t R, const slibreal_t Z) {
    bool a =
        (R < domain_rmin || R > domain_rmax ||
         Z < domain_zmin || Z > domain_zmax);

    outside_domain_flag = a;
    return (!a);
}

/**
 * Checks whether the given time value is equal to
 * or above the specified maximum time. This function
 * returns true only when the maximum time has been
 * reached, and false otherwise.
 *
 * The argument supplied to this function should be
 * given in seconds. It is then automatically converted
 * to whichever time unit the maxtime is specified in.
 *
 * t: Time to compare to maxtime (in seconds).
 * r: Radial component of particle position.
 * z: Vertical component of particle position.
 *
 * RETURNS true if t >= maxtime; false otherwise.
 */
bool ParticlePusher::MaxTimeReached(slibreal_t t, slibreal_t z) {
    switch (this->timeunit) {
        case ORBITTIMEUNIT_POLOIDAL:
            return (GetPoloidalTime(z) >= this->maxtime);
        case ORBITTIMEUNIT_SECONDS:
            return (t >= this->maxtime);
        default:
            throw ParticlePusherException("Internal: ParticlePusher::MaxTimeReached(): Unrecognized time unit: %d.", this->timeunit);
    }
}

/**
 * Evaluate a secondary orbit (i.e. an orbit that will be
 * used to calculate the Jacobian for the primary orbit)
 * 
 * p:  Particle defining initial state.
 * nd: Direction in which to nudge particle in order
 *     to evaluate this orbit.
 */
void ParticlePusher::EvaluateSecondaryOrbit(Particle *p, enum Particle::nudge_direction nd) {
    // Override timekeeping to get all necessary timesteps
    enum orbittimeunit tunit = this->timeunit;
    slibreal_t maxtime = this->maxtime;

    this->timeunit = ORBITTIMEUNIT_SECONDS;
    this->maxtime = this->integrator1->LastTime();

    p->Nudge(this->nudge_value, nd);
    RunIntegrator(this->equation, this->integrator2, p);
    // Restore particle
    p->Nudge(this->nudge_value, (nd==Particle::NUDGE_OUTWARDS?Particle::NUDGE_INWARDS:Particle::NUDGE_OUTWARDS));

    // Restore timekeeping
    this->timeunit = tunit;
    this->maxtime = maxtime;
}

/**
 * Run the timing integrator (i.e. the orbit that will be
 * used to calculate the poloidal transit time for a
 * particle orbit)
 *
 * p:  Particle defining the initial state.
 */
slibreal_t ParticlePusher::RunTimingIntegrator(Particle *p) {
    // Set temporary timing state
    this->timeunit = ORBITTIMEUNIT_POLOIDAL;

    RunIntegrator(this->timingEquation, this->timingIntegrator, p);

    // Reset original state
    this->timeunit = ORBITTIMEUNIT_SECONDS;
    this->maxtime  = this->timingIntegrator->LastTime();

    return this->maxtime;
}

/**
 * Push the given particle, i.e. calculate the
 * particle's orbit.
 *
 * p: Particle to evaluate orbit for.
 *
 * RETURNS an Orbit object representing the orbit
 * of the pushed particle.
 */
Orbit *ParticlePusher::Push(Particle *p) {
    slibreal_t poltime=0.0;

    // Run timing integrator if necessary to determine
    // the poloidal transit time
    if (this->useTimingIntegrator)
        RunTimingIntegrator(p);

    RunIntegrator(this->equation, this->integrator1, p);

    // If the chosen time unit is 'poloidal time', and
    // the particle appears to not move in the poloidal
    // plane, then this particle cannot be pushed any further
    if (restflag) {
        orbit_class_t cl = ORBIT_CLASS_STAGNATION;
        slibreal_t time = this->integrator1->LastTime();

        return retorbit->Create(time, this->integrator1, nullptr, this->timingIntegrator, this->equation, p, this->nudge_value, cl, this->forceNumericalJacobian);
    } else if (outside_domain_flag) {
        orbit_class_t cl = ORBIT_CLASS_COLLIDED;
        slibreal_t time = this->integrator1->LastTime();

        return retorbit->Create(time, this->integrator1, nullptr, this->timingIntegrator, this->equation, p, this->nudge_value, cl, this->forceNumericalJacobian);
    }

    if (this->timeunit == ORBITTIMEUNIT_POLOIDAL) {
        if (this->useTimingIntegrator)
            poltime = FindPoloidalTime(this->timingEquation, this->timingIntegrator);
        else
            poltime = FindPoloidalTime(this->equation, this->integrator1);
    }

    orbit_class_t cl1;
    if (outside_domain_flag)
        cl1 = ORBIT_CLASS_COLLIDED;
    else
        cl1 = this->equation->ClassifyOrbit(this->integrator1);

    if (this->calculateJacobianOrbit && (!this->magfield->HasMagneticFlux() || this->forceNumericalJacobian)) {
        EvaluateSecondaryOrbit(p, Particle::NUDGE_OUTWARDS);

        if (outside_domain_flag) {
            slibreal_t time = this->integrator1->LastTime();
            return retorbit->Create(time, this->integrator1, nullptr, this->timingIntegrator, this->equation, p, this->nudge_value, ORBIT_CLASS_COLLIDED, this->forceNumericalJacobian);
        }

        orbit_class_t cl2 = this->equation->ClassifyOrbit(this->integrator2);

        if (cl1 != cl2)
            EvaluateSecondaryOrbit(p, Particle::NUDGE_INWARDS);
    }

    // Choose last time step
    switch (this->timeunit) {
        case ORBITTIMEUNIT_SECONDS:
            if (this->calculateJacobianOrbit && (!this->magfield->HasMagneticFlux() || this->forceNumericalJacobian))
                return retorbit->Create(this->maxtime, this->integrator1, this->integrator2, this->timingIntegrator, this->equation, p, this->nudge_value, cl1, this->forceNumericalJacobian);
            else
                return retorbit->Create(this->maxtime, this->integrator1, nullptr, this->timingIntegrator, this->equation, p, this->nudge_value, cl1, this->forceNumericalJacobian);
        default:
            if (this->calculateJacobianOrbit && (!this->magfield->HasMagneticFlux() || this->forceNumericalJacobian))
                return retorbit->Create(poltime, this->integrator1, this->integrator2, this->timingIntegrator, this->equation, p, this->nudge_value, cl1, this->forceNumericalJacobian);
            else
                return retorbit->Create(poltime, this->integrator1, nullptr, this->timingIntegrator, this->equation, p, this->nudge_value, cl1, this->forceNumericalJacobian);
    }
}

/**
 * Initialize and take all steps with the given
 * integrator. The integrator must have already been
 * initialized with the given equation.
 */
void ParticlePusher::RunIntegrator(SOFTEquation *eqn, Integrator<6> *intg, Particle *p) {
    Vector<6> init;
    eqn->InitializeParticle(p, init);
    intg->InitialValue(init);

    this->rprev = this->rinit = GetPositionR(eqn, intg);
    this->zprev = this->zinit = GetPositionZ(eqn, intg);

    ResetDomainCheck();
    ResetPoloidalTime();

    // Advance in time
    slibreal_t T, R, Z;
    do {
        T = intg->Step();
        R = GetPositionR(eqn, intg);
        Z = GetPositionZ(eqn, intg);
    } while (!MaxTimeReached(T, Z) && InsideDomainBounds(R, Z));
}

/**
 * Determines whether 'integrator2' should be used or not.
 * The purpose of 'integrator2' is to solve for an orbit which
 * is slightly offset from the requested one, in order to
 * calculate the spatial Jacobian determinant, which is required
 * when integrating over real space.
 *
 * j: If true, runs 'integrator2' to solve for an extra orbit
 *    that can be used to calculate the spatial Jacobian.
 */
void ParticlePusher::EnableJacobianCalculation() {
    this->calculateJacobianOrbit = true;
}
