/**
 * Definition of the SOFT class.
 */
#ifndef _SOFT_H
#define _SOFT_H

#include <string>
#include <vector>

/* SOFTLIB stuff */
#include <softlib/config.h>
#include <softlib/DistributionFunction/DistributionFunction.h>
#include <softlib/MagneticField/MagneticField2D.h>

/* SOFT stuff */
#include "config.h"
#include "Orbit/ParticlePusher.h"
#include "PhaseSpace/ParticleGenerator.h"
#include "Tools/Tool.h"

struct global_settings {
    std::string distribution;
	bool include_drifts;
    std::string magnetic_field;
    unsigned int num_threads;
    std::string particle_generator;
    std::string particle_pusher;
    std::vector<std::string> tools;
};

class SOFT {
	private:
		struct global_settings *globset;
	public:
        Configuration *configuration;
        DistributionFunction *distribution;
		MagneticField2D *magfield;
		ParticleGenerator *partgen;
        //ParticlePusher *pusher;

        SOFT();
        ~SOFT();

        void Run();
        //void RunParticles();

        // Getters
        struct global_settings *GetGlobalSettings() { return this->globset; }

        // Setters
        void SetGlobalSettings(struct global_settings *g) { this->globset = g; }

        // Output routines
        typedef enum {
            MESSAGE_GENERAL,
            // Tools/ParticleGenerator
            WARNING_PG_INEFFICIENT_PARALLELIZATION,
            // Orbit/ParticlePusher
            WARNING_OPP_UNRECOGNIZED_SETTING,
            // Tools/Radiation
            WARNING_TR_NO_OUTPUT,
            // Tools/Radiation/Models/AngularDistribution
            WARNING_TRMAD_ABNORMALLY_LARGE_BETA,
            // Tools/Radiation/Models/Cone
            WARNING_TRMC_BUGGY_CONE_MODEL,
            // Tools/Radiation/Output/Green
            WARNING_TROG_IMAGE_NOT_LAST,
            // [Must always come last]
            MESSAGE_LAST
        } message_t;

        static bool VerifyMessage(const message_t);

        template<typename ... Args>
        static void PrintError(const message_t, const std::string&, Args&& ...);
        template<typename ... Args>
        static void PrintError(const std::string&, Args&& ...);
        template<typename ... Args>
        static void PrintWarning(const message_t, const std::string&, Args&& ...);
        template<typename ... Args>
        static void PrintWarning(const std::string&, Args&& ...);
        static void PrintInfo();
        template<typename ... Args>
        static void PrintInfo(const message_t, const std::string&, Args&& ...);
        template<typename ... Args>
        static void PrintInfo(const std::string&, Args&& ...);

#ifdef WITH_MPI
        template<typename ... Args>
        static void PrintMPI(const std::string&, Args&& ...);
#endif

        static const std::string PRINT_YES, PRINT_NO;

    private:
        static bool message_checklist[MESSAGE_LAST];
};

// Implementation of template members
#include "SOFT.tcc"

#endif/*_SOFT_H*/
