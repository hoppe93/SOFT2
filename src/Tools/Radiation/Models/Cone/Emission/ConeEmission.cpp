/**
 * Implementation of the 'ConeEmission' base class.
 *
 * ON THE PRESENCE OF I, Q, U AND V
 *   Note that while the parent class 'Cone' also allocates
 *   arrays for storing the four Stokes parameters, they
 *   are designed for storing the total **detected**
 *   radiation, whereas the arrays for storing Stokes
 *   parameters here are designed to store the total
 *   **emitted** radiation.
 */

#include "Tools/Radiation/Models/Cone/ConeEmission.h"

using namespace __Radiation;

/**
 * Constructor.
 * Sets up the four Stokes parameter arrays.
 *
 * det: Detector used to measure the emitted radiation.
 */
ConeEmission::ConeEmission(Detector *det) {
    this->detector = det;
    this->nwavelengths = det->GetNWavelengths();

    if (nwavelengths > 0) {
        this->wavelengths = det->GetWavelengths();
        this->I = new slibreal_t[nwavelengths];
        this->Q = new slibreal_t[nwavelengths];
        this->U = new slibreal_t[nwavelengths];
        this->V = new slibreal_t[nwavelengths];

        for (unsigned int i = 0; i < nwavelengths; i++) {
            I[i] = Q[i] = U[i] = V[i] = 0.0;
        }
    } else {
        this->I = this->Q = this->U = this->V = nullptr;
    }
}

/**
 * Destructor.
 */
ConeEmission::~ConeEmission() {
    if (this->I != nullptr)
        delete [] this->I;
    if (this->Q != nullptr)
        delete [] this->Q;
    if (this->U != nullptr)
        delete [] this->U;
    if (this->V != nullptr)
        delete [] this->V;
}
