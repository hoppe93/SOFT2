/**
 * Methods for the main SOFT loop
 */

#include <omp.h>
#include "SOFT.h"
#include "SOFTException.h"
#include "SOFTLocal.h"

/**
 * Run SOFT
 */
void SOFT::Run() {
    SOFTLocal *gs;
    bool haserror = false;
    unsigned int total_invalid = 0;

    // Initialize local SOFT object on root thread.
    gs = new SOFTLocal(this, 0, true);

    #pragma omp parallel num_threads(globset->num_threads) shared(total_invalid)
    {
        unsigned int thread_num, invalid;
        SOFTLocal *s;

        thread_num = omp_get_thread_num();
        try {
            // Clone objects on all but the root thread.
            if (thread_num == 0)
                s = gs;
            else
                s = new SOFTLocal(this, thread_num, true);

            s->Run(&invalid);
            s->Finish();

            #pragma omp atomic
            total_invalid += invalid;
        } catch (SOFTLibException& ex) {
            SOFT::PrintError(ex.whats());
            haserror = true;
        }
    }

    /* Orbits launched on the magnetic axis are typically
     * discarded. At most, there will be 'Np1 * Np2' such
     * orbits. If there are more invalid orbits, that's
     * a sign that something is weird. */
    if (total_invalid > this->partgen->GetN1()*this->partgen->GetN2())
        SOFT::PrintWarning("%u/%u orbits were considered invalid.", total_invalid, this->partgen->Size());

    if (haserror) return;

    try {
        gs->Output();
    } catch (SOFTException& ex) {
        SOFT::PrintError(ex.whats());
    }
}

/**
 * Finish the simulation.
 */
void SOFTLocal::Finish() {
    this->thandler->Finish();
}

/**
 * Generate output.
 */
void SOFTLocal::Output() {
    this->thandler->Output();
}

/**
 * Iterate over all of phase-space, generating an orbit
 * for each of the "particles" and then passing the
 * results to the available tools.
 *
 * invalid: On return, contains the number of orbits
 *          that were considered invalid in this run.
 */
void SOFTLocal::Run(unsigned int *invalid) {
    Particle *p;
    ParticlePusher *pp;
    Orbit *o;
    unsigned int inv = 0;

    p = this->partgen->AllocateParticle();
    pp = this->pusher;

    if (this->id == 0)
        Welcome();

    // Main SOFT loop
    while (!this->partgen->IsFinished()) {
        if (!this->partgen->Generate(p, this->magfield, this->distribution))
            break;
        
        o = pp->Push(p);

        /* Invalid orbit (for example, poloidal
           time is used and particle is on the
           magnetic axis) */
        if (o == nullptr) {
            inv++;
            continue;
        }

        // Pass orbit to tools (whether it was completed or not)
        this->thandler->Handle(o, p);
    }

    *invalid = inv;
}

