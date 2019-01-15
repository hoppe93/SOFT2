/**
 * Configure the 'Radiation' tool.
 */


#include <cstring>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "PhaseSpace/ParticleGenerator.h"
#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/Detector.h"

// Models
#include "Tools/Radiation/Models/AngularDistribution.h"
#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Isotropic.h"

// Output
#include "Tools/Radiation/Output/Green.h"
#include "Tools/Radiation/Output/Image.h"
#include "Tools/Radiation/Output/SoVVolume.h"
#include "Tools/Radiation/Output/Space3D.h"
#include "Tools/Radiation/Output/Spectrum.h"
#include "Tools/Radiation/Output/Topview.h"

using namespace std;
using namespace __Radiation;

int Radiation::CONFBLOCK_T_DETECTOR;
int Radiation::CONFBLOCK_T_MODEL;
int Radiation::CONFBLOCK_T_OUTPUT;

/**
 * HOW TO ADD A NEW RADIATION MODEL/OUTPUT
 *
 * 1. Add an entry to the 'radmodels/radoutputs' array below.
 *    Each entry has two properties:
 *      [0]: Tool name (configuration block name).
 *      [1]: A pointer to the function InitModel
 *           with an explicit type specification for
 *           the model Model class you're adding.
 *           For the 'Isotropic' model this would be
 *
 *             InitModel<Isotropic>
 * 2. Increase the 'RADIATION_NMODELS/RADIATION_NOUTPUTS'
 *    variable by one.
 */
template<class T>
Model *InitModel(struct global_settings*, ConfigBlock*, ConfigBlock*, Radiation*);
template<class T>
RadiationOutput *InitRadiationOutput(ConfigBlock*, ConfigBlock*, Detector*, MagneticField2D*, ParticleGenerator*);

const unsigned int RADIATION_NMODELS=3;
struct radiation_modelspec radmodels[RADIATION_NMODELS] = {
    {"isotropic", InitModel<Isotropic>},
    {"cone", InitModel<Cone>},
    {"angdist", InitModel<AngularDistribution>}
};

const unsigned int RADIATION_NOUTPUTS=6;
struct radiation_outputspec radoutputs[RADIATION_NOUTPUTS] = {
    {"green", InitRadiationOutput<Green>},
    {"image", InitRadiationOutput<Image>},
    {"sovvolume", InitRadiationOutput<SoVVolume>},
    {"space3d", InitRadiationOutput<Space3D>},
    {"spectrum", InitRadiationOutput<Spectrum>},
    {"topview", InitRadiationOutput<Topview>}
};

/**
 * Initialize a new Model object from the
 * given configuration. The model type to
 * initialize must be explicitly specified.
 *
 * conf: ConfigBlock configuring the Model object.
 * root: Root ConfigBlock, for access to ConfigBlocks
 *       for any sub-modules.
 */
template<class T>
Model *InitModel(struct global_settings *globset, ConfigBlock *conf, ConfigBlock *root, Radiation *parent) {
    T *t = new T(parent);
    t->Configure(globset, conf, root);

    return t;
}

/**
 * Initialize a new RadiationOutput object from the
 * given configuration. The output type to
 * initialize must be explicitly specified.
 *
 * conf: ConfigBlock configuring the RadiationOutput object.
 * root: Root ConfigBlock, for access to ConfigBlocks
 *       for any sub-modules.
 */
template<class T>
RadiationOutput *InitRadiationOutput(ConfigBlock *conf, ConfigBlock *root, Detector *d, MagneticField2D *m, ParticleGenerator *pgen) {
    T *t = new T(d, m, pgen);
    t->Configure(conf, root);

    return t;
}

/**
 * Set up this Radiation tool object.
 *
 * globset: Global SOFT settings.
 * conf:    Configuration block for this object.
 * root:    Root configuration block, possible containing
 *          configuration blocks for submodules used
 *          by this tool.
 */
void Radiation::Configure(
    struct global_settings *globset,
    ConfigBlock *conf, ConfigBlock *root
) {
    this->SetName(conf->GetName());

    // detector
    if (!conf->HasSetting("detector"))
        throw RadiationException("No detector specified.");
    else {
        string detname = (*conf)["detector"];
        if (!root->HasSubBlock(CONFBLOCK_T_DETECTOR, detname))
            throw RadiationException("No detector with name '%s' defined.", detname.c_str());
        
        ConfigBlock *cb = root->GetConfigBlock(CONFBLOCK_T_DETECTOR, detname);
        SetDetector(new __Radiation::Detector(cb));
    }

    // ntoroidal
    if (!conf->HasSetting("ntoroidal"))
        ntoroidal = 3500;
    else {
        Setting *s = conf->GetSetting("ntoroidal");
        if (!s->IsUnsignedInteger32())
            throw RadiationException("Invalid format of toroidal resolution 'ntoroidal'. Expected positive integer.");
        else ntoroidal = s->GetUnsignedInteger32();

        if (ntoroidal == 0)
            throw RadiationException("Invalid value assigned to 'ntoroidal': %u.", ntoroidal);
    }

    // Set up toroidal increments
    dphi    = 2.0*M_PI / ntoroidal;
    //cosdphi = cos(dphi);
    //sindphi = sin(dphi);
    cosphi = new slibreal_t[ntoroidal];
    sinphi = new slibreal_t[ntoroidal];
    for (unsigned int i = 0; i < ntoroidal; i++) {
        cosphi[i] = cos(i*dphi);
        sinphi[i] = sin(i*dphi);
    }

    // model
    if (!conf->HasSetting("model"))
        throw RadiationException("No radiation model specified.");
    else {
        string modname = (*conf)["model"];
        if (!root->HasSubBlock(CONFBLOCK_T_MODEL, modname))
            throw RadiationException("No radiation model with name '%s' defined.", modname.c_str());

        SetModel(SetupRadiationModel(globset, root->GetConfigBlock(CONFBLOCK_T_MODEL, modname), root));
    }

    // output
    if (!conf->HasSetting("output"))
        throw RadiationException("No radiation output handler(s) specified.");
    else {
        Setting *s = conf->GetSetting("output");
        if (s->GetNumberOfValues() == 0)
            throw RadiationException("No radiation output handler(s) specified.");

        vector<string> outputs = s->GetTextVector();
        this->noutput = outputs.size();

        if (this->output != nullptr)
            delete [] this->output;

        this->output = new RadiationOutput*[noutput];

        for (size_t i = 0; i < outputs.size(); i++) {
            string oname = outputs[i];
            if (!root->HasSubBlock(CONFBLOCK_T_OUTPUT, outputs[i]))
                throw RadiationException("No output module named '%s' defined.", oname.c_str());
            ConfigBlock *cb = root->GetConfigBlock(CONFBLOCK_T_OUTPUT, outputs[i]);

            this->output[i] = SetupRadiationOutput(cb, root);
            this->measuresPolarization |= this->output[i]->MeasuresPolarization();
        }
    }

    // torthreshold
    if (conf->HasSetting("torthreshold")) {
        Setting *s = conf->GetSetting("torthreshold");
        if (!s->IsScalar())
            throw RadiationException("Invalid value assigned to 'torthreshold'. Expected real value.");

        this->torthreshold = s->GetScalar();

        if (this->torthreshold < 0 || this->torthreshold > 1)
            throw RadiationException("Invalid value assigned to 'torthreshold'. Must be between 0 and 1.");
    }

    // torquad
    if (conf->HasSetting("torquad")) {
        if ((*conf)["torquad"] == "maximize")
            this->quadrature = Radiation::QUADRATURE_FINDSOV;
        else if ((*conf)["torquad"] == "trapz")
            this->quadrature = Radiation::QUADRATURE_TRAPEZOIDAL;
        else
            throw RadiationException("Invalid number of values assigned to option 'quadrature'.");
    } else
        this->quadrature = Radiation::QUADRATURE_FINDSOV;

    // Initialize flags used in improved trapz
    if (this->quadrature == Radiation::QUADRATURE_FINDSOV) {
        this->torflags = new char[ntoroidal];
        memset(this->torflags, 0, ntoroidal);
    }

    // wall opacity
    if (conf->HasSetting("wall_opacity")) {
        Setting *s = conf->GetSetting("wall_opacity");
        if (s->GetNumberOfValues() != 1)
            throw RadiationException("Invalid value assigned to 'wall_opacity'.");

        if (s->GetString() == "opaque")
            this->wall_opacity = WALL_OPACITY_OPAQUE;
        else if (s->GetString() == "semi")
            this->wall_opacity = WALL_OPACITY_SEMI_TRANSPARENT;
        else if (s->GetString() == "transparent")
            this->wall_opacity = WALL_OPACITY_TRANSPARENT;
        else
            throw RadiationException("Unrecognized setting for 'wall_opacity': %s.", s->GetString().c_str());
    } else
        this->wall_opacity = WALL_OPACITY_SEMI_TRANSPARENT;
}

/**
 * Prepare a configuration object for reading
 * settings for this tool.
 *
 * conf: Configuration object to prepare.
 */
void Radiation::PrepareConfiguration(Configuration *conf) {
    CONFBLOCK_T_DETECTOR = conf->RegisterBlockType("@Detector");
    CONFBLOCK_T_MODEL    = conf->RegisterBlockType("@RadiationModel");
    CONFBLOCK_T_OUTPUT   = conf->RegisterBlockType("@RadiationOutput");
}

/**
 * Configure the radiation model to use.
 *
 * conf: Configuration block for the radiation model.
 * root: Root ConfigBlock (for access to ConfigBlocks of
 *       possible sub-modules).
 */
Model *Radiation::SetupRadiationModel(struct global_settings *globset, ConfigBlock *conf, ConfigBlock *root) {
    string stype = conf->GetSecondaryType();

    for (unsigned int i = 0; i < RADIATION_NMODELS; i++) {
        if (radmodels[i].name == stype) {
            return radmodels[i].init(globset, conf, root, this);
        }
    }

    throw RadiationException("No radiation model of type '%s' available in this version of SOFT.", stype.c_str());
}

/**
 * Configure all radiation outputs to use.
 *
 * conf: Configuration block for the radiation output.
 * root: Root ConfigBlock (for access to ConfigBlocks of
 *       possible sub-modules).
 */
RadiationOutput *Radiation::SetupRadiationOutput(ConfigBlock *conf, ConfigBlock *root) {
    string stype = conf->GetSecondaryType();

    for (unsigned int i = 0; i < RADIATION_NOUTPUTS; i++) {
        if (radoutputs[i].name == stype) {
            return radoutputs[i].init(conf, root, detector, magfield, pgen);
        }
    }

    throw RadiationException("No radiation output module of type '%s' available in this version of SOFT.", stype.c_str());
}

