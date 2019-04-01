#ifndef _PARTICLE_GENERATOR_H
#define _PARTICLE_GENERATOR_H

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/ProgressTracker.h>

class ParticleGenerator;

#include "PhaseSpace/Particle.h"
#include "SOFT.h"

class ParticleGenerator {
	private:
		slibreal_t r0, r1, dr,
				   p10, p11, dp1,
				   p20, p21, dp2,
                   **rhoeff, rhomax, rhomin;
        slibreal_t *rgrid=nullptr, *p1grid=nullptr, *p2grid=nullptr;
		slibreal_t charge, mass;
		unsigned int ir, i1, i2,
					 nr, n1, n2;
		int mom1type, mom2type;				// Momentum types

		int specified_position=Particle::POSITION_PARTICLE;

		slibreal_t drift_shift_tolerance=1e-4;	// Desired tolerance when computing orbit drift shift

		/* Flags */
		bool include_drifts=false;			// Account for drifts when calculating radial position of particles
		bool finished=false;				// true if all particles have been generated. false otherwise.
        bool print_progress=false;          // Print phase space progress continuously

        ProgressTracker *progress;

        enum MPI_Distribute_Mode {
            MPI_DISTMODE_ALL,
            MPI_DISTMODE_RADIUS,
            MPI_DISTMODE_MOMENTUM1,
            MPI_DISTMODE_MOMENTUM2
        };

        enum MPI_Distribute_Mode mpi_distribute_mode=MPI_DISTMODE_ALL;
	public:
		ParticleGenerator(MagneticField2D*, ConfigBlock*, struct global_settings*);
        ~ParticleGenerator();

		Particle *AllocateParticle();
		slibreal_t CalculateRadialOrbitDriftShift(MagneticField2D*, const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t);
        slibreal_t CalculateVerticalOrbitDriftShift(MagneticField2D*, const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t);
		void _Calculate_Xpol(MagneticField2D*, bool, const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t, const slibreal_t, slibreal_t&, slibreal_t&);

		bool Generate(Particle*, MagneticField2D*, DistributionFunction *f=nullptr);
        void GenerateCoordinateGrids(int mpi_distribute_parameter=Particle::COORDINATE_);
        void GenerateRhoeffTable(MagneticField2D*);
        void GetRadialCoordinate(ConfigBlock*, const std::string&, slibreal_t, slibreal_t);
		void InitializeParticle(
            Particle*, DistributionFunction*, MagneticField2D*,
            const slibreal_t, const slibreal_t, const slibreal_t,
            const unsigned int, const unsigned int, const unsigned int
        );
		bool IsFinished();
        void SetRadialCoordinate(MagneticField2D*, ConfigBlock*);

        // Getters
        unsigned int GetNr() { return this->nr; }
        unsigned int GetN1() { return this->n1; }
        unsigned int GetN2() { return this->n2; }
        unsigned int Size()  { return this->nr*this->n1*this->n2; }

        int GetP1Type() { return mom1type; }
        int GetP2Type() { return mom2type; }

        slibreal_t *GetRGrid()  { return this->rgrid; }
        slibreal_t *GetP1Grid() { return this->p1grid; }
        slibreal_t *GetP2Grid() { return this->p2grid; }

        slibreal_t GetDriftShiftTolerance() { return this->drift_shift_tolerance; }

        enum MPI_Distribute_Mode GetMPIDistributionMode() { return this->mpi_distribute_mode; }
};

class ParticleGeneratorException : public SOFTException {
    public:
        template<typename ... Args>
        ParticleGeneratorException(const std::string &msg, Args&& ... args)
            : SOFTException(msg, std::forward<Args>(args) ...) {
            AddModule("Particle generator");
        }
};

#endif/*_PARTICLE_GENERATOR_H*/
