/**
 * PARTICLE GENERATOR
 * Implementation of the ParticleGenerator class, that takes a grid definition
 * as input and generates a succession of 'Particle' objects that should be
 * run by SOFT.
 *
 * In difference to many other objects, this object has no default settings to
 * prevent the user from using several incompatible momentum-space coordinates
 * when defining the grid.
 *
 * USAGE:
 *    
 *    // Globally
 *    ParticleGenerator *pg =
 *        new ParticleGenerator(magnetic_field, conf);
 *
 *    ...
 *
 *    // On each thread...
 *    Particle *p;
 *    p = pg->AllocateParticle();
 *
 *    ...
 *
 *    // In loop, on each thread
 *    // Generate one particle
 *    if (pg->Generate(p)) {
 *        // Run particle
 *        ...
 *    } else {
 *        // Finished
 *    }
 */

#include <string>
#include <functional>

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>

#include "config.h"
#include "constants.h"
#include "PhaseSpace/Particle.h"
#include "PhaseSpace/ParticleGenerator.h"
#include "SOFT.h"
#include "units.h"

#ifdef WITH_MPI
#   include <mpi.h>
#endif

using namespace std;

/**
 * Constructor.
 *
 * mf:   Pointer to SOFT magnetic field object.
 * conf: Configuration object defining the grid.
 */
const int pg_ncoordinates=7;
const string pg_coordinate_names[pg_ncoordinates] = { "gamma", "p", "ppar", "pperp", "thetap", "ithetap", "xi" };
const int pg_coordinate_types[pg_ncoordinates] = {
	Particle::COORDINATE_GAMMA,
	Particle::COORDINATE_P,
	Particle::COORDINATE_PPAR,
	Particle::COORDINATE_PPERP,
	Particle::COORDINATE_THETAP,
	Particle::COORDINATE_ITHETAP,
	Particle::COORDINATE_XI
};
ParticleGenerator::ParticleGenerator(MagneticField2D *mf, ConfigBlock *conf, struct global_settings *glob) {
#	define PG_MAX_COORDINATES 2
	int coordinates[PG_MAX_COORDINATES], ncoords=0, i;
	vector<slibreal_t> coordvals[2], radius;
	Setting *s;

    ir = i1 = i2 = 0;

	// Copy interesting global settings
	this->include_drifts = glob->include_drifts;

	// Gather momentum-space coordinates
	for (i = 0; i < pg_ncoordinates; i++) {
		if (conf->HasSetting(pg_coordinate_names[i])) {
			if (ncoords == PG_MAX_COORDINATES)
				throw ParticleGeneratorException("Too many momentum-space coordinates provided to particle generator. Give exactly two coordinates.");

			s = conf->GetSetting(pg_coordinate_names[i]);
			if (!s->IsNumericVector())
				throw ParticleGeneratorException("Invalid specification of coordinate '%s'. Expected numeric vector with three (3) elements.", pg_coordinate_names[i].c_str());

			coordinates[ncoords] = pg_coordinate_types[i];
			coordvals[ncoords++] = s->GetNumericVector();
		}
	}
	
	if (ncoords != PG_MAX_COORDINATES)
		throw ParticleGeneratorException("Too few momentum-space coordinates provided to particle generator. Give exactly two coordinates.");
	
	// Make sure valid syntax was used for coordinate vectors
	for (i = 0; i < 2; i++) {
		if (coordvals[i].size() != 3)
			throw ParticleGeneratorException(
				"Coordinate %s: Expected exactly arguments to the coordinate. Syntax: %s=[start],[end],[number of points].",
				pg_coordinate_names[coordinates[i]].c_str(), pg_coordinate_names[coordinates[i]].c_str()
			);
	}

	// Set momentum coordinates
	this->mom1type = coordinates[0];
	this->mom2type = coordinates[1];
	this->p10 = coordvals[0][0];
	this->p11 = coordvals[0][1];
	this->n1  = (unsigned int)coordvals[0][2];
	this->p20 = coordvals[1][0];
	this->p21 = coordvals[1][1];
	this->n2  = (unsigned int)coordvals[1][2];

    if (this->n1 > 1)
        this->dp1 = (this->p11-this->p10) / ((slibreal_t)(this->n1-1));
    else if (this->n1 < 1)
        throw ParticleGeneratorException(
            "Invalid number of points specified for velocity parameter '%s': %u.\n",
            pg_coordinate_names[this->mom1type], this->n1
        );
    else
        this->dp1 = 0.0;

    if (this->n2 > 1)
        this->dp2 = (this->p21-this->p20) / ((slibreal_t)(this->n2-1));
    else if (this->n2 < 1)
        throw ParticleGeneratorException(
            "Invalid number of points specified for velocity parameter '%s': %u.\n",
            pg_coordinate_names[this->mom2type], this->n2
        );
    else
        this->dp2 = 0.0;

    // Verify that coordinates are given valid values.
    // (these trow exceptions if not).
    Particle::VerifyCoordinateSpecification(this->p10, this->p11, this->n1, this->mom1type);
    Particle::VerifyCoordinateSpecification(this->p20, this->p21, this->n2, this->mom2type);
	
	// Get radial coordinate
    SetRadialCoordinate(mf, conf);

#ifdef WITH_MPI
    if (conf->HasSetting("mpi_distribute_mode")) {
        s = conf->GetSetting("mpi_distribute_mode");

        if (s->GetNumberOfValues() != 1)
            throw ParticleGeneratorException(
                "Too many parameters specified for 'mpi_distribute_mode'. Expected exactly one parameter."
            );

        string str = s->GetString();
        if (str == "radius" || str == "a" || str == "r" || str == "rho")
            this->mpi_distribute_mode = MPI_DISTMODE_RADIUS;
        else if (str == "1")
            this->mpi_distribute_mode = MPI_DISTMODE_MOMENTUM1;
        else if (str == "2")
            this->mpi_distribute_mode = MPI_DISTMODE_MOMENTUM2;
        else {
            for (int i = 0; i < pg_ncoordinates; i++) {
                if (str == pg_coordinate_names[i]) {
                    if (pg_coordinate_types[i] == this->mom1type)
                        this->mpi_distribute_mode = MPI_DISTMODE_MOMENTUM1;
                    else if (pg_coordinate_types[i] == this->mom2type)
                        this->mpi_distribute_mode = MPI_DISTMODE_MOMENTUM2;
                    else
                        throw ParticleGeneratorException(
                            "Unable to distribute parameter '%s' over MPI. The parameter is not part of the phase space.",
                            str.c_str()
                        );

                    break;
                }
            }

            if (i == pg_ncoordinates)
                throw ParticleGeneratorException(
                    "Unrecognized parameter '%s' specified to be distributed over MPI.",
                    str.c_str()
                );
        }
    }
#endif

    GenerateCoordinateGrids(this->mpi_distribute_mode);

    // Get mass and charge (if set)
    if (conf->HasSetting("charge")) {
        s = conf->GetSetting("charge");
        this->charge = s->GetScalar() * CHARGE_UNIT;
    } else this->charge = -1.0 * CHARGE_UNIT;

    if (conf->HasSetting("mass")) {
        s = conf->GetSetting("mass");
        if (s->GetScalar() <= 0.0)
            throw ParticleGeneratorException("Invalid mass assigned to particle. Mass must be positive.");
        this->mass = s->GetScalar() * MASS_UNIT;
    } else this->mass = 1.0 * MASS_UNIT;

    // What position is specified? (particle or GC?)
    if (conf->HasSetting("position")) {
        s = conf->GetSetting("position");
        if (s->GetString() == "guiding-center" ||
            s->GetString() == "gc")
            this->specified_position = Particle::POSITION_GUIDINGCENTER;
        else if (s->GetString() == "particle")
            this->specified_position = Particle::POSITION_PARTICLE;
        else
            throw ParticleGeneratorException("Unrecognized choice for specified positition: '"+s->GetString()+"'.");
    } else
        this->specified_position = Particle::POSITION_GUIDINGCENTER;

    // Other settings
    if (conf->HasSetting("driftshifttol")) {
        s = conf->GetSetting("driftshifttol");
        if (s->IsScalar() && s->GetScalar() < 1.0 && s->GetScalar() > REAL_EPSILON)
            this->drift_shift_tolerance = s->GetScalar();
        else
            throw ParticleGeneratorException("Invalid value assigned to 'driftshifttol'.");
    }

    if (conf->HasSetting("progress")) {
        size_t nprint_progress = 10;
        ProgressTracker::ProgressType ptype = ProgressTracker::PROGRESS_LINES;
        const bool ESTIMATE_PROGRESS = true;

        s = conf->GetSetting("progress");

        if (s->IsBool(0))
            print_progress = s->GetBool(0);
        else if (s->IsUnsignedInteger32(0)) {
            print_progress = true;
            nprint_progress = s->GetUnsignedInteger32(0);

            if (nprint_progress == 0)
                throw ParticleGeneratorException("Invalid value assigned to 'progress'. Must be greater than zero.");
        } else
            throw ParticleGeneratorException("Invalid value assigned to 'progress'.");

        if (s->GetNumberOfValues() == 2) {
            if (s->GetString(1) == "lines")
                ptype = ProgressTracker::PROGRESS_LINES;
            else
                throw ParticleGeneratorException("Unrecognized progress type assigned to 'progress': '%s'.", s->GetString(1).c_str());
        } else if (s->GetNumberOfValues() > 2)
            throw ParticleGeneratorException("Too many values assigned to 'progress'. Expected one (1) or two (2) values.");

        if (print_progress) {
            size_t total = ((size_t)nr) * ((size_t)n1) * ((size_t)n2);

#ifdef COLOR_TERMINAL
            progress = new ProgressTracker(total, nprint_progress, ptype, true, ESTIMATE_PROGRESS);
#else
            progress = new ProgressTracker(total, nprint_progress, ptype, false, ESTIMATE_PROGRESS);
#endif
        }
    } else
        print_progress = false;

    // Generate table of effective magnetic axis location
    // (if drifts are enabled)
    if (this->include_drifts)
        GenerateRhoeffTable(mf);
}
ParticleGenerator::~ParticleGenerator() {
    delete [] this->p2grid;
    delete [] this->p1grid;
    delete [] this->rgrid;
}

/**
 * Generate the individual coordinate grids of the phase space.
 */
void ParticleGenerator::GenerateCoordinateGrids(
    enum MPI_Distribute_Mode
#ifdef WITH_MPI
    param
#endif
) {
    unsigned int i;

    this->end_ir = this->nr;
    this->end_i1 = this->n1;
    this->end_i2 = this->n2;

    // TODO
#ifdef WITH_MPI
    int nprocesses, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (param == MPI_DISTMODE_RADIUS) {
        unsigned int dr = nr / nprocesses;
        if (mpi_rank >= (nr % nprocesses))
            dr++;

        this->ir     = mpi_rank*dr;
        this->end_ir = ir + dr;
    } else if (param == MPI_DISTMODE_MOMENTUM1) {
        unsigned int d1 = n1 / nprocesses;
        if (mpi_rank >= (n1 % nprocesses))
            d1++;

        this->i1     = mpi_rank*d1;
        this->end_i1 = i1 + d1;
    } else if (param == MPI_DISTMODE_MOMENTUM2) {
        unsigned int d2 = n2 / nprocesses;
        if (mpi_rank >= (n2 % nprocesses))
            d2++;

        this->i2     = mpi_rank*d2;
        this->end_i2 = i2 + d2;
    } else
        throw ParticleGeneratorException(
            "GenerateCoordinateGrids(): Unrecognized MPI distribution mode: %d.", param
        );
#endif

    rgrid  = new slibreal_t[this->nr];
    p1grid = new slibreal_t[this->n1];
    p2grid = new slibreal_t[this->n2];

    for (i = 0; i < this->nr; i++)
        rgrid[i]  = this->r0  + i * this->dr;
    for (i = 0; i < this->n1; i++)
        p1grid[i] = this->p10 + i * this->dp1;
    for (i = 0; i < this->n2; i++)
        p2grid[i] = this->p20 + i * this->dp2;
}

/**
 * Sets the radial coordinate of the particle generator using
 * the given configuration block. The radial coordinate can be
 * specified in three ways: (1) major radius, (2) minor radius
 * or (3) normalized minor radius. Internally, SOFT only works
 * in (3).
 */
void ParticleGenerator::SetRadialCoordinate(MagneticField2D *mf, ConfigBlock *conf) {
    this->rhomin = mf->GetMagneticAxisR();
    this->rhomax = mf->GetMaxRadius();

    if (conf->HasSetting("a"))
        GetRadialCoordinate(conf, "a", this->rhomin, this->rhomax - this->rhomin);
    else if (conf->HasSetting("r"))
        GetRadialCoordinate(conf, "r", this->rhomin, 1.0);
    else if (conf->HasSetting("rho"))
        GetRadialCoordinate(conf, "rho", 0.0, 1.0);
    else
        throw ParticleGeneratorException("No radial coordinate specified.");

    if (this->r0 < this->rhomin)
        throw ParticleGeneratorException("Invalid value assigned to inner initial radius. Must be on or to the right of the magnetic axis.");
    else if (this->r1 > this->rhomax)
        SOFT::PrintWarning("Outer initial radius is outside domain.");
        //throw ParticleGeneratorException("Invalid value assigned to outer initial radius. Must be inside device/separatrix.");
}

/**
 * Reads the radial coordinate with name 'name' from the
 * given ConfigBlock and sets the radial parameters
 * r0, r1 and nr. They are set to
 *
 *    r0 = c + scale * setting[0]
 *    r1 = c + scale * setting[1]
 *    nr = setting[2]
 *
 * conf:  ConfigBlock object containing the setting.
 * name:  Name of radial coordinate to read.
 * c:     Constant by which to shift radius.
 * scale: Factor by which to scale input radius.
 */
void ParticleGenerator::GetRadialCoordinate(
    ConfigBlock *conf, const string &name,
    slibreal_t c, slibreal_t scale
) {
    Setting *s;

	if (conf->HasSetting(name)) {			// Normalized minor radius
		s = conf->GetSetting(name);
        if (s->GetNumberOfValues() == 1) {
            if (!s->IsUnsignedInteger32())
                throw ParticleGeneratorException("Invalid specification of coordinate '%s'. Expected '[[scalar,]scalar,]uint32'.", name.c_str());

            this->r0 = this->rhomin;
            this->r1 = this->rhomax;
            this->nr = s->GetUnsignedInteger32();
        } else if (s->GetNumberOfValues() == 2) {
            if (!s->IsScalar(0) || !s->IsUnsignedInteger32(1))
                throw ParticleGeneratorException("Invalid specification of coordinate '%s'. Expected '[[scalar,]scalar,]uint32'.", name.c_str());

            this->r0 = this->rhomin;
            this->r1 = c + s->GetScalar(0) * scale;
            this->nr = s->GetUnsignedInteger32(1);
        } else if (s->GetNumberOfValues() == 3) {
            if (!s->IsScalar(0) || !s->IsScalar(1) || !s->IsUnsignedInteger32(2))
                throw ParticleGeneratorException("Invalid specification of coordinate '%s'. Expected '[[scalar,]scalar,]uint32'.", name.c_str());

            this->r0 = c + s->GetScalar(0) * scale;
            this->r1 = c + s->GetScalar(1) * scale;
            this->nr = s->GetUnsignedInteger32(2);
        } else
            throw ParticleGeneratorException("Invalid specification of coordinate '%s'. Expected '[[scalar,]scalar,]uint32'.", name.c_str());
	} else
        throw ParticleGeneratorException("Radial coordinate '%s' not specified.", name.c_str());

    // Verify order
    if (this->r0 > this->r1)
        throw ParticleGeneratorException("The initial radius may not be larger than the final radius.");

    // Set radial step
    if (this->nr < 1)
        throw ParticleGeneratorException("Invalid number of radial points specified: %u.", this->nr);
    else if (this->nr == 1) {
        if (this->r0 != this->r1)
            throw ParticleGeneratorException("nr = 1, but initial and final radii are different.");

        this->dr = 0.0;
    } else
        this->dr = (this->r1-this->r0) / (this->nr-1);
}

/**
 * Generate a table of values for the effective magnetic
 * axis as a function of the two momentum coordinates.
 */
void ParticleGenerator::GenerateRhoeffTable(
    MagneticField2D *magfield
) {
    slibreal_t p1, p2, ppar, pperp;
    unsigned int i, j;

    rhoeff    = new slibreal_t*[n1];
    rhoeff[0] = new slibreal_t[n1*n2];

    for (i = 0; i < n1; i++) {
        p1 = p10 + i*dp1;

        if (i > 0)
            rhoeff[i] = rhoeff[i-1] + n2;

        for (j = 0; j < n2; j++) {
            p2 = p20 + j*dp2;
            Particle::ToPP(p1, p2, this->mom1type, this->mom2type, &ppar, &pperp);

            rhoeff[i][j] = this->rhomin + CalculateRadialOrbitDriftShift(
                magfield, this->mass, this->charge, ppar, pperp
            );
        }
    }
}

/**
 * Allocates and initializes a Particle object
 * that can later be passed to 'Generate()' to
 * generate a new particle.
 *
 * RETURNS a Particle object initialized with
 *    the most basic properties (charge, drift
 *    shift, mass and position type).
 */
Particle *ParticleGenerator::AllocateParticle() {
	Particle *p = new Particle();

	p->SetCharge(this->charge);
	p->SetDriftShift(0.0);
	p->SetMass(this->mass);
	p->SetPositionType(this->specified_position);

	return p;
}

/**
 * Initialize the given particle with the given
 * coordinate values.
 *
 * part: Particle object to initialize.
 * f:    Distribution function from which to draw the particle.
 * mf:   Magnetic field in which particle lives.
 * rho:  Radial location of particle.
 * p1:   Momentum coordinate 1.
 * p2:   Momentum coordinate 2.
 * ir:   Radial coordinate index.
 * i1:   Momentum coordinate 1 index.
 * i2:   Momentum coordinate 2 index.
 */
void ParticleGenerator::InitializeParticle(
	Particle *part, DistributionFunction *f,
    MagneticField2D *mf,
    const slibreal_t rho, const slibreal_t p1,
    const slibreal_t p2, const unsigned int ir,
    const unsigned int i1, const unsigned int i2
) {
	slibreal_t d=0.0, _dp1=(dp1==0?1:dp1), _dp2=(dp2==0?1:dp2), z0;
    part->SetIndices(ir, i1, i2);
	part->InitializeMomentum(mom1type, mom2type, p1, p2, _dp1, _dp2);

	if (include_drifts)
        d = this->rhoeff[i1][i2] - this->rhomin;

    z0 = this->CalculateVerticalOrbitDriftShift(
        mf, this->mass, this->charge, part->GetPpar(), part->GetPperp(), rho
    );

    // Evaluate distribution function
	if (f != nullptr)
		part->SetF(f->Eval(rho, part->GetMomentum(), part->GetXi(), d));
	else
		part->SetF(1.0);

    if (dr == 0)
        part->InitializePosition(this->specified_position, rho, z0, 1.0, d);
    else
        part->InitializePosition(this->specified_position, rho, z0, dr, d);
}

/**
 * Generates the next particle in line to be
 * simulated.
 *
 * part: The particle object in which to store
 *       the generated particle.
 * f:    Distribution function from which to draw
 *       the particle.
 *
 * RETURNS 'true' if a new particle was generated.
 *    'false' if there are no more particles in the
 *    queue. Note that if 'false' is returned,
 *    the 'part' object has NOT been populated
 *    with any new data.
 */
bool ParticleGenerator::Generate(Particle *part, MagneticField2D *mf, DistributionFunction *f) {
	slibreal_t r, p1, p2;
	bool success = true;
    unsigned int lir=ir, li1=i1, li2=i2;

	#pragma omp critical (ParticleGenerator_Generate)
	{
		if (this->finished) success = false;
		else {
            r  = this->rgrid[ir];
            p1 = this->p1grid[i1];
            p2 = this->p2grid[i2];

            lir = ir;
            li1 = i1;
            li2 = i2;

			ir++;

			if (ir >= end_ir) {
				ir = 0;
				i1++;
				
				if (i1 >= end_i1) {
					i1 = 0;
					i2++;

					if (i2 >= end_i2) {
						this->finished = true;
						i2 = 0;
					}
				}
			}
		}
	}

    if (print_progress) {
        size_t indx = ((size_t)lir) + ((size_t)nr)*(((size_t)li1) + ((size_t)n1)*((size_t)li2));
        progress->PrintProgress(indx+1);
    }

	if (success) {
        if (!include_drifts || r >= this->rhoeff[li1][li2])
            this->InitializeParticle(part, f, mf, r, p1, p2, lir, li1, li2);
        else
            return Generate(part, mf, f);
    }

	return success;
}

/**
 * Checks if the particle generator is done
 * generating all particles in the queue.
 *
 * RETURNS true if there are no more particles
 * to generate. false otherwise.
 */
bool ParticleGenerator::IsFinished() { return this->finished; }

/**
 * Calculates the expected shift of the
 * effective magnetic axis due to orbit drifts
 * in the radial direction. This is done using
 * a basic line search algorithm, to minimize
 * the poloidal guiding-center speed.
 *
 * ppar:  Parallel momentum of particle.
 * pperp: Perpendicular momentum of particle.
 *
 * RETURNS the amount by which the magnetic
 * axis is shifted in the radial direction
 * (including sign) due to drifts.
 */
slibreal_t ParticleGenerator::CalculateRadialOrbitDriftShift(
    MagneticField2D *magfield, const slibreal_t m, const slibreal_t q,
    const slibreal_t ppar, const slibreal_t pperp
) {
	const slibreal_t invphi = 0.5*(sqrt(5.0)-1.0),	// Inverse golden ratio
        invphi2 = invphi*invphi;
	slibreal_t a = magfield->GetMinRadius(),
			   b = magfield->GetMaxRadius(),
			   c, d, h,
			   rmaj = magfield->GetMagneticAxisR(),
               zaxis = magfield->GetMagneticAxisZ(),
			   Xpol_c = 0, Xpol_d = 0;
	
    auto xpol = [this, &magfield, &m, &q, &ppar, &pperp, &zaxis](const slibreal_t x) {
        slibreal_t Xdotr, Xdotz;
        _Calculate_Xpol(magfield, true, m, q, ppar, pperp, x, zaxis, Xdotr, Xdotz);
        return hypot(Xdotr, Xdotz);
    };

    h = b-a;
	c = a + h*invphi2;
	d = a + h*invphi;
	Xpol_c = xpol(c);
	Xpol_d = xpol(d);

    unsigned int n = ceil(log(this->drift_shift_tolerance/h) / log(invphi)), i;

    for (i = 0; i < n; i++) {
		if (Xpol_c < Xpol_d) {
			b = d;
			d = c;
			Xpol_d = Xpol_c;
            h *= invphi;
			c = a + h*invphi2;
			Xpol_c = xpol(c);
		} else {
			a = c;
			c = d;
			Xpol_c = Xpol_d;
            h *= invphi;
			d = a + h*invphi;
			Xpol_d = xpol(d);
		}
	}

    if (Xpol_c < Xpol_d)
        return ((a+d)*0.5 - rmaj);
    else
        return ((b+c)*0.5 - rmaj);
}

/**
 * Calculates the expected shift of the
 * effective magnetic axis due to orbit drifts
 * in the radial direction. This is done using
 * a basic line search algorithm, to minimize
 * the poloidal guiding-center speed.
 *
 * ppar:  Parallel momentum of particle.
 * pperp: Perpendicular momentum of particle.
 *
 * RETURNS the amount by which the magnetic
 * axis is shifted in the radial direction
 * (including sign) due to drifts.
 */
slibreal_t ParticleGenerator::CalculateVerticalOrbitDriftShift(
    MagneticField2D *magfield, const slibreal_t m, const slibreal_t q,
    const slibreal_t ppar, const slibreal_t pperp, const slibreal_t r
) {
	const slibreal_t invphi = 0.5*(sqrt(5.0)-1.0),	// Inverse golden ratio
        invphi2 = invphi*invphi;
	slibreal_t a=0, b=0,
			   c, d, h,
			   Xpol_c = 0, Xpol_d = 0;
	
    if (r == magfield->GetMagneticAxisR())
        return magfield->GetMagneticAxisZ();

    auto xpol = [this, &magfield, &m, &q, &ppar, &pperp, &r](const slibreal_t z) {
        slibreal_t Xdotr, Xdotz;
        _Calculate_Xpol(magfield, this->include_drifts, m, q, ppar, pperp, r, z, Xdotr, Xdotz);

        return fabs(Xdotr);
    };

    // Bracket the interval
    a = magfield->GetMagneticAxisZ();
    b = a + SQRT_REAL_EPSILON;

    Xpol_c = xpol(a);
    Xpol_d = xpol(b);
    c = xpol(a-SQRT_REAL_EPSILON);

    d = Xpol_d - Xpol_c;        // f(x+h) - f(x)
    slibreal_t db = 0.05 * (r-magfield->GetMagneticAxisR()) * (d>0?(+1):(-1));

    if (d > 0) {
        b = a - db;
        Xpol_d = xpol(b);
    }

    while ((Xpol_d=xpol(b)) < Xpol_c)
        b -= db;
    
    // Carry out the golden-section search
    if (a > b) {
        c = a;
        a = b;
        b = c;
    }

    h = b-a;

    if (h <= this->drift_shift_tolerance)
        return 0.5*(a+b);

    c = a + h*invphi2;
    d = a + h*invphi;
    Xpol_c = xpol(c);
    Xpol_d = xpol(d);

    unsigned int n = ceil(log(this->drift_shift_tolerance/h) / log(invphi)), i;


    for (i = 0; i < n; i++) {
		if (Xpol_c < Xpol_d) {
			b = d;
			d = c;
			Xpol_d = Xpol_c;
            h *= invphi;
			c = a + h*invphi2;
			Xpol_c = xpol(c);
		} else {
			a = c;
			c = d;
			Xpol_c = Xpol_d;
            h *= invphi;
			d = a + h*invphi;
			Xpol_d = xpol(d);
		}
	}

    if (Xpol_c < Xpol_d)
        return ((a+d)*0.5);
    else
        return ((b+c)*0.5);
}

/**
 * INTERNAL ROUTINE
 * Calculates the guiding-center drift speed if the
 * particle has parallel and perpendicular momentum
 * ppar and pperp respectively, in the radial point r.
 *
 * magfield: Magnetic field to evaluate Xpol in.
 * wdrifts:  Take drifts into account.
 * m:        Particle mass.
 * q:        Particle charge.
 * ppar:     Parallel momentum of particle.
 * pperp:    Perpendicular momentum of particle.
 * r:        Radial coordinate of particle.
 * z:        Vertical coordinate of particle.
 *
 * RETURNS
 * Xdotr:    Radial component of guiding-center velocity.
 * Xdotz:    Vertical component of guiding-center velocity.
 */
void ParticleGenerator::_Calculate_Xpol(
    MagneticField2D *magfield, bool wdrifts,
    const slibreal_t m, const slibreal_t q,
    const slibreal_t ppar, const slibreal_t pperp, const slibreal_t r,
    const slibreal_t z, slibreal_t &Xdotr, slibreal_t &Xdotz
) {
	struct magnetic_field_data mfd = magfield->EvalDerivatives(r, 0, z);
	slibreal_t Beffpar;
    const slibreal_t
        c  = LIGHTSPEED,
        mu = m*c*c*pperp*pperp / (2.0*mfd.Babs);

    if (wdrifts) {
        Vector<3> Beff, B(mfd.B), bhat(B/mfd.Babs), gradB(mfd.gradB), curlB(mfd.curlB);
        Vector<3> bhatXgradB = Vector<3>::Cross(bhat, gradB);
        
        Beff = B + (m*c*ppar/(q*mfd.Babs)) * (curlB + bhatXgradB);
        Beffpar = Beff.Dot(bhat);

        Xdotr = (c*ppar*Beff[0] + mu / q * bhatXgradB[0]) / Beffpar;
        Xdotz = (c*ppar*Beff[2] + mu / q * bhatXgradB[2]) / Beffpar;
    } else {
        Vector<3> B(mfd.B), bhat(B / mfd.Babs);
        Xdotr = bhat[0];
        Xdotz = bhat[1];
    }
}

