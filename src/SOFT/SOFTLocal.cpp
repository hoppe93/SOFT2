/**
 * Implementation of the 'local' SOFT object.
 *
 * This SOFT object is designed to run asynchronously with other
 * objects of the same type. It is created by a parent 'SOFT' class
 * which has been configured appropriately.
 */

#include <softlib/MagneticField/MagneticField2D.h>
#include "SOFT.h"
#include "SOFTLocal.h"

/**
 * Constructor.
 *
 * parent: Parent SOFT object.
 * clone:  If true, clones all relevant objects in the parent SOFT
 *         object. If false, all pointers to relevant objects are
 *         instead copied. The latter option should only be used for
 *         one object. Default: true.
 */
SOFTLocal::SOFTLocal(SOFT *parent, unsigned int id, bool clone) {
    this->id = id;
    this->soft = parent;
    this->partgen = parent->partgen;

    /* We clone these objects rather than creating
     * them from scratch because:
     *
     *   distribution: The 'MinClone()' available for distributions allows
     *                 the memory usage to be optimized. Also, we won't
     *                 need to load a numerical distribution from disk
     *                 several times.
     *   magfield:     Loading a numeric magnetic field from the hard drive
     *                 can be unnecessarily slow to do in parallel.
     */
    if (clone) {
        this->magfield = parent->magfield->Clone();
        this->distribution = parent->distribution->MinClone();
    } else {
        this->magfield = parent->magfield;
        this->distribution = parent->distribution;
    }

    this->pusher = InitPusher(this->magfield, parent->GetGlobalSettings(), parent->configuration);
    this->thandler = InitTools(parent, parent->GetGlobalSettings(), parent->configuration, this->partgen, this->pusher, this->magfield);
}
