/**
 * Radiation :: Output :: Space3D
 * 
 * Handle output from a RadiationModel.
 */

#include <omp.h>
#include <softlib/config.h>
#include "Tools/Radiation/Output/Space3D.h"

using namespace __Radiation;
using namespace std;

/**
 * Handle output from the given model.
 *
 * det:   Detector observing the radiation.
 * model: Radiation model to handle output from.
 * rp:    Particle producing the radiation.
 */
void Space3D::Handle(Detector *__UNUSED__(det), Model *model, RadiationParticle *rp) {
    Vector<3> IJK = rp->GetPosition() - point0;
    Vector<3> dif = point1 - point0;

    long long int I = (long long int)((IJK[0]/dif[0]) * pixelsX);
    long long int J = (long long int)((IJK[1]/dif[1]) * pixelsY);
    long long int K = (long long int)((IJK[2]/dif[2]) * pixelsZ);

    // Is pixel within image?
    if (I < 0 || I >= pixelsX ||
        J < 0 || J >= pixelsY ||
        K < 0 || K >= pixelsZ)
        return;

    size_t index = (I*pixelsY + J)*pixelsZ + K;

    slibreal_t f = rp->GetF();
    slibreal_t diffel = rp->GetDifferentialElement();

    slibreal_t p = model->GetPower() * diffel*f;

    #pragma omp atomic
    s3dimage[index] += p;
}

