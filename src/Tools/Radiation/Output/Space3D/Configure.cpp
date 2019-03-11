/**
 * Radiation :: Output :: Space3D
 *
 * Configuration of the 'Space3D' radiation output module.
 */

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Tools/Radiation/Output/Space3D.h"

using namespace __Radiation;
using namespace std;

/**
 * Configure the 'Space3D' radiation output module.
 * 
 * conf: Configuration block for the module.
 * root: Root configuration block, providing access
 *       to the ConfigBlock's of possible sub-modules (unused).
 */
void Space3D::Configure(ConfigBlock *conf, ConfigBlock *__UNUSED__(root)) {
    this->SetName(conf->GetName());

    // output
    if (!conf->HasSetting("output"))
        throw Space3DException("No S3D output name has been specifed. Set the option 'output'.");
    else
        this->output = (*conf)["output"];
    
    // pixels
    if (!conf->HasSetting("pixels"))
        throw Space3DException("The number of pixels of S3D was not specified.");
    else {
        Setting *s = conf->GetSetting("pixels");
        if (!s->IsUnsignedInteger64())
            throw Space3DException("Invalid value assigned to parameter 'pixels'.");

        this->pixelsX = s->GetUnsignedInteger64();
        this->pixelsY = this->pixelsZ = this->pixelsX;

        if (this->pixelsX == 0)
            throw Space3DException("Invalid value assigned to parameter 'pixels'.");
    }

    // point0
    GetPoint(conf, "point0", point0);
    GetPoint(conf, "point1", point1);

    this->imagesize = pixelsX*pixelsY*pixelsZ;
}

/**
 * Set the 'point' parameter with name 'name' from
 * the configuration block 'conf'. Assign the value
 * to 'p'.
 *
 * conf: Configuration block to load value from.
 * name: Name of point to read.
 * p:    Vector representing point to assign.
 */
void Space3D::GetPoint(ConfigBlock *conf, const string &name, Vector<3> &p) {
    if (!conf->HasSetting(name))
        throw Space3DException("The parameter '%s' was not configured.", name.c_str());

    Setting *s = conf->GetSetting(name);

    if (!s->IsNumericVector(3))
        throw Space3DException(
            "Invalid value assigned to parameter '%s'. Expected 3 scalars.",
            name.c_str()
        );

    p[0] = s->GetScalar(0);
    p[1] = s->GetScalar(1);
    p[2] = s->GetScalar(2);
}

/**
 * Initialize the 3D image.
 */
void Space3D::Initialize() {
    #pragma omp critical (Space3D_Initialize)
    {
        if (s3dimage == nullptr) {
            s3dimage = new slibreal_t[imagesize];

            for (size_t i = 0; i < imagesize; i++)
                s3dimage[i] = 0.0;
        }
    }
}

