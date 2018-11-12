/**
 * Radiation :: Output :: Topview
 * 
 * Handle output from a RadiationModel.
 */

#include <softlib/config.h>
#include "Tools/Radiation/Output/Topview.h"

using namespace __Radiation;

/**
 * Calculate the pixel in which the given
 * particle is observed.
 *
 * det: Detector to consider.
 * rp:  Particle emitting radiation.
 * i:   Contains horizontal pixel index on return.
 * j:   Contains vertical pixel index on return.
 */
void Topview::GetTopviewPixel(RadiationParticle *rp, unsigned int *i, unsigned int *j) {
    Vector<3> &X = rp->GetPosition();

    *i = lround(X[0] / max_radius * npixels * 0.5) + npixels/2;
    *j = lround(X[1] / max_radius * npixels * 0.5) + npixels/2;
}

/**
 * Handle output from the given model.
 *
 * det:   Detector observing the radiation.
 * model: Radiation model to handle output from.
 * rp:    Particle producing the radiation.
 */
void Topview::Handle(Detector *__UNUSED__(det), Model *model, RadiationParticle *rp) {
    unsigned int i, j;
    slibreal_t diffel, f;

    GetTopviewPixel(rp, &i, &j);
    if (i >= npixels || j >= npixels)
        return;

    f = rp->GetF();
    diffel = rp->GetDifferentialElement();

    this->topview[i*npixels + j] += model->GetPower() * diffel * f;
}

