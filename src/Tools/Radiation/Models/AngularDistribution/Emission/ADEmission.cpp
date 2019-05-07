/**
 * Implementation of the base class 'ADEmission'.
 */

#include <softlib/config.h>
#include "SOFT.h"
#include "Tools/Radiation/Models/AngularDistribution/ADEmission.h"

using namespace __Radiation;

/**
 * Constructor.
 *
 * det:        Detector object defining the detector properties.
 */
ADEmission::ADEmission(Detector *det, MagneticField2D *mf) {
    this->detector = det;
    this->magfield = mf;
    this->nwavelengths = det->GetNWavelengths();

    if (this->nwavelengths > 0) {
        this->wavelengths = det->GetWavelengths();
        this->I = new slibreal_t[this->nwavelengths];
        this->Q = new slibreal_t[this->nwavelengths];
        this->U = new slibreal_t[this->nwavelengths];
        this->V = new slibreal_t[this->nwavelengths];
    } else {
        this->wavelengths = nullptr;
        this->I = this->Q = this->U = this->V = nullptr;
    }
}

/**
 * Destructor.
 */
ADEmission::~ADEmission() {
    if (this->I != nullptr)
        delete [] this->I;
    if (this->Q != nullptr)
        delete [] this->Q;
    if (this->U != nullptr)
        delete [] this->U;
    if (this->V != nullptr)
        delete [] this->V;
}

