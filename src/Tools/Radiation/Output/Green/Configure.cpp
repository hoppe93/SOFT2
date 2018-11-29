/**
 * Radiation :: Output :: Green
 *
 * Configuration of the 'Green' radiation output module.
 */

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "SOFT.h"
#include "Tools/Radiation/Output/Green.h"

using namespace std;
using namespace __Radiation;

/**
 * Allocate memory for the image.
 * 
 * pixelsset: True if the number of pixels has been set.
 *            False otherwise.
 */
void Green::AllocateGreen(bool pixelsset) {
    int i;

    if (this->storeStokesParameters)
        this->fsize = 4;
    else
        this->fsize = 1;

    // Set up Green's function "factors"
    // (the "factors" array is an array consisting of the factors
    // which we should multiply each index with to retrieve the
    // corresponding element in the Green's function; thus it
    // should be built by multiplying the dimensions of the
    // Green's function cumulatively)
    this->factors = new size_t[this->nformat];
    for (i = 0; i < this->nformat; i++)
        this->factors[i] = 1;

    size_t s;
    for (i = 0; i < this->nformat; i++) {
        switch (this->format[i]) {
            case '1': s = this->n1; this->hasP1 = true; this->i1 = i; break;
            case '2': s = this->n2; this->hasP2 = true; this->i2 = i; break;
            case 'i': s = this->subnrowpixels; this->hasI = true; this->ii = i; break;
            case 'j': s = this->subncolpixels; this->hasJ = true; this->ij = i; break;
            case 'r': s = this->nr; this->hasR = true; this->ir = i; break;
            case 'w': s = this->nw; this->hasW = true; this->iw = i; break;
            default:
                throw GreenException("Invalid character in Green's function format string: %c.", this->format[i]);
        }

        this->fsize *= s;

        if (i > 0)
            this->factors[i-1] = s;
    }

    this->fsizeWithoutStokes = this->fsize / 4;

    // Actually construct factors array
    for (int j = this->nformat-2; j >= 0; j--)
        this->factors[j] *= this->factors[j+1];

    if ((this->hasI || this->hasJ) && !pixelsset)
        throw GreenException("Green's function contains pixels, but the number of pixels was not set.");

    this->function = new slibreal_t[this->fsize];
    this->containsAllMomentumSpaceParameters = (this->hasP1 && this->hasP2 && this->hasR);

    // Initialize function
    for (size_t i = 0; i < this->fsize; i++)
        this->function[i] = 0;
}

/**
 * Configure the 'Green' radiation output module.
 * 
 * conf: Configuration block for the module.
 * root: Root configuration block, providing access
 *       to the ConfigBlock's of possible sub-modules (unused).
 */
void Green::Configure(ConfigBlock *conf, ConfigBlock *__UNUSED__(root)) {
    Setting *s;
    bool pixelsset = false;

    this->SetName(conf->GetName());

    // format
    if (!conf->HasSetting("format"))
        throw GreenException("No Green's function format has been specified.");
    else {
        this->format = (*conf)["format"];
        ValidateFormat(this->format);

        this->nformat = this->format.length();
    }

    // output
    if (!conf->HasSetting("output"))
        throw GreenException("No Green's function output name has been specifed. Set the option 'output'.");
    else
        this->output = (*conf)["output"];

    // pixels
    if (conf->HasSetting("pixels")) {
        s = conf->GetSetting("pixels");
        if (s->GetNumberOfValues() == 1)
            this->nrowpixels = this->ncolpixels = s->GetUnsignedInteger32();
        else if (s->GetNumberOfValues() == 2) {
            this->nrowpixels = s->GetUnsignedInteger32(0);
            this->ncolpixels = s->GetUnsignedInteger32(1);
        } else
            throw GreenException("Invalid specification of 'pixels'.");

        pixelsset = true;
    }

    // stokesparams
    if (conf->HasSetting("stokesparams")) {
        s = conf->GetSetting("stokesparams");
        if (!s->IsBool())
            throw GreenException("Invalid value assigned to parameter 'stokesparams'. Expected boolean.");

        this->storeStokesParameters = s->GetBool();
    } else
        this->storeStokesParameters = false;

    // suboffseti
    if (conf->HasSetting("suboffseti")) {
        s = conf->GetSetting("suboffseti");
        if (!s->IsUnsignedInteger32())
            throw GreenException("Invalid value assigned to parameter 'suboffseti'. Expected non-negative integer.");

        this->suboffseti = s->GetUnsignedInteger32();
    } else
        this->suboffseti = 0;

    // suboffsetj
    if (conf->HasSetting("suboffsetj")) {
        s = conf->GetSetting("suboffsetj");
        if (!s->IsUnsignedInteger32())
            throw GreenException("Invalid value assigned to parameter 'suboffsetj'. Expected non-negative integer.");

        this->suboffsetj = s->GetUnsignedInteger32();
    } else
        this->suboffsetj = 0;

    // subpixels
    if (conf->HasSetting("subpixels")) {
        s = conf->GetSetting("subpixels");
        if (s->GetNumberOfValues() == 1)
            this->subnrowpixels = this->subncolpixels = s->GetUnsignedInteger32();
        else if (s->GetNumberOfValues() == 2) {
            this->subnrowpixels = s->GetUnsignedInteger32(0);
            this->subncolpixels = s->GetUnsignedInteger32(1);
        } else
            throw GreenException("Invalid specification of 'subpixels'.");
    } else if (pixelsset) {
        this->subnrowpixels = this->nrowpixels;
        this->subncolpixels = this->ncolpixels;
    }

    // with_jacobian
    if (conf->HasSetting("with_jacobian")) {
        s = conf->GetSetting("with_jacobian");
        if (!s->IsBool())
            throw GreenException("Invalid value assigned to 'with_jacobian'. Expected 'yes' or 'no'.");

        this->withJacobian = s->GetBool();
    } else
        this->withJacobian = true;

    if (pixelsset)
        ValidateSubPixels();

    AllocateGreen(pixelsset);
}

/**
 * Validate the given string to make sure that
 * it has a valid format. Also warns about
 * potentially bad combinations of parameters.
 *
 * fmt: String to validate.
 */
void Green::ValidateFormat(const string &fmt) {
    #define NTOUCHED 128
    size_t n = fmt.length(), i;
    bool touched[NTOUCHED] = {false};
    int c;

    for (i = 0; i < n; i++) {
        c = fmt[i];
        if (c < 0 || c > NTOUCHED-1)
            throw GreenException("Invalid character in Green's function format at position %zu.", i+1);

        if (touched[c])
            throw GreenException("Parameter '%c' used multiple times in Green's function format.", c);

        switch (c) {
            case '1': case '2': case 'i':
            case 'j': case 'r': case 'w':
                touched[c] = true;
                break;
            default:
                throw GreenException("Unrecognized character in Green's function format: %c.", c);
        }
    }

    if (touched['i']) {
        if (!touched['j'])
            throw GreenException("Only column pixel index specified to Green's function. Both row and column pixels must be specified simultaneously.");
        else if (fmt.substr(n-2, 2) != "ij" && fmt.substr(n-2, 2) != "ji")
            SOFT::PrintWarning(SOFT::WARNING_TROG_IMAGE_NOT_LAST, "It is recommended to put image (ij) at the end of the Green's function.");
    } else if (touched['j'])
        throw GreenException("Only column pixel index specified to Green's function. Both row and column pixels must be specified simultaneously.");
}

/**
 * Make sure that the 'subpixels' have been
 * configured properly.
 */
void Green::ValidateSubPixels() {
    int
        ti = this->subnrowpixels+this->suboffseti,
        tj = this->subncolpixels+this->suboffsetj;

    if (this->suboffseti >= this->nrowpixels)
        throw GreenException("The row pixel offset is outside image. %u >= %u.", this->suboffseti, this->nrowpixels);
    else if (this->suboffsetj >= this->ncolpixels)
        throw GreenException("The column pixel offset is outside image. %u >= %u.", this->suboffsetj, this->ncolpixels);
    else if (this->subnrowpixels > this->nrowpixels)
        throw GreenException("The number of subset row pixels are more than the total number of row pixels. %u > %u.", this->subnrowpixels, this->nrowpixels);
    else if (this->subncolpixels > this->ncolpixels)
        throw GreenException("The number of subset column pixels are more than the total number of column pixels. %u > %u.", this->subncolpixels, this->ncolpixels);
    else if (ti > this->nrowpixels)
        throw GreenException("The subset image is larger than the total image. Outermost (row) subpixel = %u > %u.", ti, this->nrowpixels);
    else if (tj > this->ncolpixels)
        throw GreenException("The subset image is larger than the total image. Outermost (col) subpixel = %u > %u.", tj, this->ncolpixels);
}

