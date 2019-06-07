/**
 * Radiation :: Output :: Green
 *
 * Configuration of the 'Green' radiation output module.
 */

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "config.h"
#include "SOFT.h"
#include "Tools/Radiation/Output/Green.h"

using namespace std;
using namespace __Radiation;

const string Green::DEFAULT_QUANTITIES[] = {
	RadiationOutput::PARAM1,
	RadiationOutput::PARAM1NAME,
	RadiationOutput::PARAM2,
	RadiationOutput::PARAM2NAME,
	RadiationOutput::R,
	RadiationOutput::RO_DOMAIN
};
template<typename T, unsigned int sz>
unsigned int __def_size(T(&)[sz]) { return sz; }
const unsigned int Green::NDEFAULT_QUANTITIES = __def_size(Green::DEFAULT_QUANTITIES);

/**
 * Allocate memory for the image.
 */
void Green::PrepareAllocateGreen() {
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

    this->ndimensions = this->nformat;
    this->dimensions = new sfilesize_t[this->nformat];

    unsigned int NR=this->nr, N1=this->n1, N2=this->n2;
#ifdef WITH_MPI
    if (this->mpi_output_mode == MPI_Output_Mode::CHUNKED) {
        NR = this->end_r - this->start_r;
        N1 = this->end_1 - this->start_1;
        N2 = this->end_2 - this->start_2;
    }
#endif

    // Update phase-space grids in case they are
    // limited by MPI "chunked" mode
    slibreal_t *gr = new slibreal_t[NR];
    slibreal_t *g1 = new slibreal_t[N1];
    slibreal_t *g2 = new slibreal_t[N2];
    
    for (unsigned int i = 0; i < NR; i++)
        gr[i] = this->rgrid[i + this->start_r];
    for (unsigned int i = 0; i < N1; i++)
        g1[i] = this->p1grid[i + this->start_1];
    for (unsigned int i = 0; i < N2; i++)
        g2[i] = this->p2grid[i + this->start_2];
    
    this->rgrid  = gr;
    this->p1grid = g1;
    this->p2grid = g2;

    size_t s;
    for (i = 0; i < this->nformat; i++) {
        switch (this->format[i]) {
            case '1': s = N1; this->hasP1 = true; this->i1 = i; break;
            case '2': s = N2; this->hasP2 = true; this->i2 = i; break;
            case 'i': s = this->subnrowpixels; this->hasI = true; this->ii = i; break;
            case 'j': s = this->subncolpixels; this->hasJ = true; this->ij = i; break;
            case 'r': s = NR; this->hasR = true; this->ir = i; break;
            case 'w': s = this->nw; this->hasW = true; this->iw = i; break;
            default:
                throw GreenException("Invalid character in Green's function format string: %c.", this->format[i]);
        }

        this->fsize *= s;
        this->dimensions[i] = s;

        if (i > 0)
            this->factors[i-1] = s;
    }

    this->fsizeWithoutStokes = this->fsize / 4;

    // Actually construct factors array
    for (int j = this->nformat-2; j >= 0; j--)
        this->factors[j] *= this->factors[j+1];

    if ((this->hasI || this->hasJ) && !this->pixelsset)
        throw GreenException("Green's function contains pixels, but the number of pixels was not set.");

#ifdef WITH_MPI
    if (this->mpi_output_mode == MPI_Output_Mode::CHUNKED) {
        if (this->mpi_distribute_mode == ParticleGenerator::MPI_DISTMODE_RADIUS && !this->hasR)
            throw GreenException("The MPI Distribute Mode is set to 'radius', but the radial parameter is not part of the Green's function.");
        else if (this->mpi_distribute_mode == ParticleGenerator::MPI_DISTMODE_MOMENTUM1 && !this->hasP1)
            throw GreenException("The MPI Distribute Mode is set to '1', but the first momentum parameter is not part of the Green's function.");
        else if (this->mpi_distribute_mode == ParticleGenerator::MPI_DISTMODE_MOMENTUM2 && !this->hasP2)
            throw GreenException("The MPI Distribute Mode is set to '2', but the second momentum parameter is not part of the Green's function.");
    }
#endif

    this->containsAllMomentumSpaceParameters = (this->hasP1 && this->hasP2 && this->hasR);
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

    this->pixelsset = false;
    this->SetName(conf->GetName());

	// common
	if (conf->HasSetting("common"))
		this->ConfigureCommonQuantities(DEFAULT_QUANTITIES, NDEFAULT_QUANTITIES, conf->GetSetting("common")->GetTextVector());
	else
		this->ConfigureCommonQuantities(DEFAULT_QUANTITIES, NDEFAULT_QUANTITIES);

    // f_as_linear_array
    if (conf->HasSetting("f_as_linear_array")) {
        s = conf->GetSetting("f_as_linear_array");
        if (!s->IsBool())
            throw GreenException("Unrecognized value assigned to 'f_as_linear_array'. Expected boolean value.");

        this->storeFAsLinearArray = s->GetBool();
    } else
        this->storeFAsLinearArray = false;

    // format
    if (!conf->HasSetting("format"))
        throw GreenException("No Green's function format has been specified.");
    else {
        this->format = (*conf)["format"];
        ValidateFormat(this->format);

        this->nformat = this->format.length();
    }

    // MPI (output) mode
    if (conf->HasSetting("mpi_mode")) {
        s = conf->GetSetting("mpi_mode");
        if (s->GetString() == "chunked")
            this->mpi_output_mode = MPI_Output_Mode::CHUNKED;
        else if (s->GetString() == "contiguous")
            this->mpi_output_mode = MPI_Output_Mode::CONTIGUOUS;
        else
            throw GreenException("Invalid value assigned to parameter 'mpi_mode'. Expected 'chunked' or 'contiguous'.");
    } else
        this->mpi_output_mode = MPI_Output_Mode::CONTIGUOUS;

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

        this->pixelsset = true;
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
    } else if (this->pixelsset) {
        this->subnrowpixels = this->nrowpixels;
        this->subncolpixels = this->ncolpixels;
    }

    // with_f
    if (conf->HasSetting("with_f")) {
        s = conf->GetSetting("with_f");
        if (!s->IsBool())
            throw GreenException("Invalid value assigned to 'with_f'. Expected 'yes' or 'no'.");

        this->weighWithDistribution = s->GetBool();
    } else
        this->weighWithDistribution = false;

    // with_jacobian
    if (conf->HasSetting("with_jacobian")) {
        s = conf->GetSetting("with_jacobian");
        if (!s->IsBool())
            throw GreenException("Invalid value assigned to 'with_jacobian'. Expected 'yes' or 'no'.");

        this->withJacobian = s->GetBool();
    } else
        this->withJacobian = true;

    if (this->pixelsset)
        ValidateSubPixels();

    PrepareAllocateGreen();
}

/**
 * Initialize this output module (after configuration).
 */
void Green::Initialize() {
    try {
        this->function = new slibreal_t[this->fsize];
    } catch (bad_alloc& ba) {
        throw GreenException("Failed to allocate memory for Green's function: %s.", ba.what());
    }

    // Initialize function
    for (size_t i = 0; i < this->fsize; i++)
        this->function[i] = 0;
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

    if (touched[(int)'i']) {
        if (!touched[(int)'j'])
            throw GreenException("Only column pixel index specified to Green's function. Both row and column pixels must be specified simultaneously.");
        else if (fmt.substr(n-2, 2) != "ij" && fmt.substr(n-2, 2) != "ji")
            SOFT::PrintWarning(SOFT::WARNING_TROG_IMAGE_NOT_LAST, "It is recommended to put image (ij) at the end of the Green's function.");
    } else if (touched[(int)'j'])
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

