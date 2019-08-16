/**
 * Configuration of the 'Cone' radiation model.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "SOFT.h"
#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Cone/ConeBremsstrahlungEmission.h"
#include "Tools/Radiation/Models/Cone/ConeBremsstrahlungScreenedEmission.h"
#include "Tools/Radiation/Models/Cone/ConeSynchrotronEmission.h"
#include "Tools/Radiation/Models/Cone/ConeUnitEmission.h"
#include "Tools/Radiation/Models/Cone/Projection/ConeProjection.h"
#include "Tools/Radiation/Models/Cone/Projection/Original.h"
#include "Tools/Radiation/Models/Cone/Projection/Reverse.h"


using namespace std;
using namespace __Radiation;

/**
 * Configure the cone model.
 *
 * conf: ConfigBlock for the model.
 * root: Root ConfigBlock. Provided to allow access
 *       to configuration of possible sub-modules.
 */
void Cone::Configure(struct global_settings *__UNUSED__(globset), ConfigBlock *conf, ConfigBlock *__UNUSED__(root)) {
    // edgecheck
    if (conf->HasSetting("edgecheck")) {
        Setting *s = conf->GetSetting("edgecheck");
        if (!s->IsBool())
            throw ConeException("Invalid value assigned to parameter 'edgecheck'. Expected boolean.");

        this->edgeCheck = s->GetBool();
    } else
        this->edgeCheck = false;

    // emission
    if (!conf->HasSetting("emission"))
        throw ConeException("No emission type specified.");
    ConfigureEmission((*conf)["emission"], conf);

    // projection
    if (conf->HasSetting("projection")) {
        if ((*conf)["projection"] == "original") {
            this->projection = new ConeProjectionOriginal(this->parent->detector);
            SOFT::PrintWarning(SOFT::WARNING_TRMC_BUGGY_CONE_MODEL, "A buggy cone model is being used. Please, be careful.");
        } else if ((*conf)["projection"] == "reverse")
            this->projection = new ConeProjectionReverse(this->parent->detector);
        else
            throw ConeException("Invalid value assigned to parameter 'projection'.");
    } else
        this->projection = new ConeProjectionReverse(this->parent->detector);
}

/**
 * Configure the emission module.
 *
 * emname: Name of the emission module to use.
 */
void Cone::ConfigureEmission(const string& emname, ConfigBlock *conf) {
    if (emname == "bremsstrahlung") {
        if (this->parent->MeasuresPolarization())
            throw ConeBremsstrahlungException("The bremsstrahlung radiation model does not support polarization measurements.");

        //default values
        slibreal_t *Z = new slibreal_t[1];
        Z[1] = 1.0;
        slibreal_t *density = new slibreal_t[1];
        density[1] = 1.0;
        unsigned int nspecies = 1;
      
    	// Reading Z's and number-desities from pi-file
        if (conf->HasSetting("Z")) {
            Setting *s = conf->GetSetting("Z");
            if (!s->IsNumericVector())
                throw ConeException("Invalid value assigned to paramter 'Z'. Expected real vector. If no Z-values are given, Z=1 by defualt.");
            vector Z_vec = s->GetNumericVector();
            nspecies = Z_vec.size(); 
            Z = new slibreal_t [nspecies];
            
            for(unsigned int i = 0; i < nspecies; i++)
                Z[i] = Z_vec[i];
        }

        if (conf->HasSetting("n")) {
            Setting *s2 = conf->GetSetting("n");
            if (!s2->IsNumericVector())
                throw ConeException("Invalid value assigned to paramter 'n'. Expected real vector. If no n-values are given, n=1 by default.");
            vector n_vec = s2->GetNumericVector();
            if(n_vec.size() != nspecies)
                throw ConeException("Number of n-values do not match number of Z-values. If no n-values are given, n=1 by defualt.");
            density = new slibreal_t [nspecies];

            for(unsigned int i = 0; i < nspecies; i++)
                density[i] = n_vec[i];
        }

        this->emission = new ConeBremsstrahlungEmission(this->parent->detector, this->parent->magfield, nspecies, Z, density); //New arguments
    } else if (emname == "bremsstrahlung_screened") {
    
        if (!conf->HasSetting("Z"))
            throw ConeException("No Z-value(s)");
        if (!conf->HasSetting("Z0"))
            throw ConeException("No Z0-value(s)");
        if (!conf->HasSetting("n"))
            throw ConeException("No number-density value(s)");
        Setting *s = conf->GetSetting("Z"),
            *s2 = conf->GetSetting("Z0"),
            *s3 = conf->GetSetting("n");

        if(!s->IsNumericVector())
            throw ConeException("Invalid value assigned to paramer Z, expected real vector");
        if(!s2->IsNumericVector())
            throw ConeException("Invalid value assigned to paramer Z0, expected real vector");
        if(!s3->IsNumericVector())
            throw ConeException("Invalid value assigned to paramer n, expected real vector");

        vector Z_vec = s->GetNumericVector(),
            Z0_vec = s2->GetNumericVector(),
            n_vec = s3->GetNumericVector();

        const unsigned int nspecies = Z_vec.size();
        if (Z0_vec.size() != nspecies)
            throw ConeException("Number of Z0-values do not match number of Z-values");
        if (n_vec.size() != nspecies)
            throw ConeException("Number of n-values do not match number of Z-values");
        
        slibreal_t *Z = new slibreal_t [nspecies],
            *Z0 = new slibreal_t [nspecies],
            *density = new slibreal_t [nspecies];

        for(unsigned int i = 0; i < nspecies; i++){ 
            Z[i] = Z_vec[i];
            Z0[i] = Z0_vec[i];
            density[i] = n_vec[i];
        }
        this->emission = new ConeBremsstrahlungScreenedEmission(this->parent->detector, this->parent->magfield, nspecies, Z, Z0, density);
    }
else if (emname == "synchrotron") {
        this->emission = new ConeSynchrotronEmission(this->parent->detector, this->parent->magfield, this->parent->MeasuresPolarization());
    } else if (emname == "unit") {
        if (this->parent->MeasuresPolarization())
            throw ConeUnitException("The 'Unit' emission model does not support polarization measurements.");

        this->emission = new ConeUnitEmission(this->parent->detector, this->parent->magfield);
    } else
        throw ConeException("Unrecognized emission model requested: '%s'.", emname.c_str());

    this->emissionName = emname;
}

/**
 * Returns a string describing this model.
 */
const string Cone::GetDescription() const {
#ifdef COLOR_TERMINAL
    string model = Model::CONE_COLOR+"Cone\x1B[0m  with ";
    if (this->emissionName == "bremsstrahlung")
        model += Model::BREMSSTRAHLUNG_COLOR+"bremsstrahlung\x1B[0m";
    else if (this->emissionName == "synchrotron")
        model += Model::SYNCHROTRON_COLOR+"synchrotron\x1B[0m";
    else
        model += Model::UNIT_EMISSION_COLOR+this->emissionName+"\x1B[0m";

    return model;
#else
    return "Cone  with "+this->emissionName;
#endif
}

