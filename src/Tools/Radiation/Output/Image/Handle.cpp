/**
 * Radiation :: Output :: Image
 * 
 * Handle output from a RadiationModel.
 */

#include <softlib/config.h>
#include "Tools/Radiation/Output/Image.h"

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
void Image::GetImagePixel(Detector *det, RadiationParticle *rp, int &i, int &j) {
    return Image::GetImagePixel(det, rp, this->nrowpixels, this->ncolpixels, i, j);
}
void Image::GetImagePixel(Detector *det, RadiationParticle *rp, int nrowpixels, int ncolpixels, int &i, int &j) {
    Vector<3> &rcp = rp->GetRCP(), &n = det->GetDirection();
    Vector<3> ltilde;
    slibreal_t ii, jj, nDotRcp, lfov;

    nDotRcp = n.Dot(rcp);
    ltilde = rcp - n*nDotRcp;
    ii = ltilde.Dot(det->GetEHat1());
    jj = ltilde.Dot(det->GetEHat2());

    lfov = sqrt(2.0) * nDotRcp * det->GetTanVisionAngleFOV();
    i = lround(nrowpixels*(ii+0.5*lfov) / lfov);
    j = lround(ncolpixels*(jj+0.5*lfov) / lfov);
}

/**
 * Handle output from the given model.
 *
 * det:   Detector observing the radiation.
 * model: Radiation model to handle output from.
 * rp:    Particle producing the radiation.
 */
void Image::Handle(Detector *det, Model *model, RadiationParticle *rp) {
    int i, j;
    slibreal_t diffel, f;

    GetImagePixel(det, rp, i, j);
    if (i >= nrowpixels || i < 0 || j >= ncolpixels || j < 0)
        return;

    diffel = rp->GetDifferentialElement();
    f = rp->GetF();

    slibreal_t v = model->GetPower() * diffel * f;
    this->image[i*this->ncolpixels + j] += v;

    if (this->MeasuresPolarization()) {
        this->imageQ[i*this->ncolpixels + j] += model->GetPowerQ() * diffel * f;
        this->imageU[i*this->ncolpixels + j] += model->GetPowerU() * diffel * f;
        this->imageV[i*this->ncolpixels + j] += model->GetPowerV() * diffel * f;
    }
}

