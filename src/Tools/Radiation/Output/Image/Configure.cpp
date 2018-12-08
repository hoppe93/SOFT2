/**
 * Radiation :: Output :: Image
 *
 * Configuration of the 'Image' radiation output module.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Tools/Radiation/Output/Image.h"

using namespace __Radiation;

/**
 * Allocate memory for the image.
 */
void Image::AllocateImage() {
    this->ntotpixels = this->nrowpixels*this->ncolpixels;
    this->image = new slibreal_t[this->ntotpixels];

    // Initialize
    for (int i = 0; i < this->ntotpixels; i++)
        this->image[i] = 0.0;
}

/**
 * Configure the 'Image' radiation output module.
 * 
 * conf: Configuration block for the module.
 * root: Root configuration block, providing access
 *       to the ConfigBlock's of possible sub-modules (unused).
 */
void Image::Configure(ConfigBlock *conf, ConfigBlock *__UNUSED__(root)) {
    Setting *s;

    this->SetName(conf->GetName());

    // output
    if (!conf->HasSetting("output"))
        throw ImageException("No image output name has been specifed. Set the option 'output'.");
    else
        this->output = (*conf)["output"];

    // pixels
    if (!conf->HasSetting("pixels"))
        throw ImageException("The number of pixels in the image was not set.");
    else {
        s = conf->GetSetting("pixels");
        if (s->GetNumberOfValues() == 1)
            this->nrowpixels = this->ncolpixels = s->GetUnsignedInteger32();
        else if (s->GetNumberOfValues() == 2) {
            this->nrowpixels = s->GetUnsignedInteger32(0);
            this->ncolpixels = s->GetUnsignedInteger32(1);
        } else
            throw ImageException("Invalid specification of 'pixels'.");
    }
}

/**
 * Initialize this output module (after configuration,
 * but before being used).
 */
void Image::Initialize() {
    AllocateImage();
}

