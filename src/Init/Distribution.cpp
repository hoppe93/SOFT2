/**
 * Initialization of the distribution function.
 */

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/DistributionFunction/DistributionFunction.h>
#include <softlib/DistributionFunction/AnalyticalAvalanche.h>
#include <softlib/DistributionFunction/ExponentialPitch.h>
#include <softlib/DistributionFunction/CODEDistributionFunction.h>
#include <softlib/DistributionFunction/LUKEDistributionFunction.h>
#include <softlib/DistributionFunction/LinearRadialProfile.h>
#include <softlib/DistributionFunction/PowerRadialProfile.h>
#include <softlib/DistributionFunction/RadialProfile.h>
#include <softlib/DistributionFunction/SOFTDistributionFunction.h>
#include <softlib/DistributionFunction/UniformRadialProfile.h>
#include <softlib/DistributionFunction/UnitDistributionFunction.h>
#include <softlib/MagneticField/MagneticField2D.h>

#include "Init/InitDistribution.h"
#include "SOFTException.h"

using namespace std;

/**
 * General initialization routine.
 * Initializes the distribution function according to the
 * given specification.
 *
 * conf: Configuration object specifying how to initialize
 *       the distribution function.
 */
DistributionFunction *InitDistributionFunction(MagneticField2D *magfield, ConfigBlock *conf, ConfigBlock *root) {
    DistributionFunction *df;

    string type = conf->GetSecondaryType();
    if (type == "avalanche")
        df = InitAvalancheDistribution(magfield, conf, root);
    else if (type == "code")
        df = InitCODEDistribution(magfield, conf, root);
    else if (type == "luke")
        df = InitLUKEDistribution(magfield, conf);
    else if (type == "norse")
        df = InitNORSEDistribution(magfield, conf, root);
    else if (type == "numerical")
        df = InitNumericalDistribution(magfield, conf);
    else if (type == "pitch")
        df = InitExponentialPitchDistribution(magfield, conf, root);
    else if (type == "unit")
        df = InitUnitDistributionFunction(magfield, conf, root);
    else
        throw SOFTException("Distribution function '%s': type: Unrecognized distribution function type: '%s'.", conf->GetName().c_str(), type.c_str());

    return df;
}

/**
 * Initialize distribution function according to the analytical
 * avalanche formula.
 */
RadialDistributionFunction *InitAvalancheDistribution(MagneticField2D *magfield, ConfigBlock *conf, ConfigBlock *root) {
    Setting *set;
    slibreal_t EHat, lnLambda, Zeff;
    RadialProfile *radprof;

    EHat     = InitAvalancheDistribution_param(conf, "EHat");
    lnLambda = InitAvalancheDistribution_param(conf, "lnLambda");
    Zeff     = InitAvalancheDistribution_param(conf, "Zeff");

    if (conf->HasSetting("radprof")) {
        set = conf->GetSetting("radprof");
        radprof = InitRadialProfile(magfield, set, root, conf->GetName());
    } else
        radprof = new UniformRadialProfile();

    // Combine...
    return new RadialDistributionFunction(
        radprof,
        new AnalyticalAvalanche(
            EHat, lnLambda, Zeff
        )
    );
}
/**
 * Routine for initializing a single parameter of the
 * analytical avalanche distribution function.
 *
 * conf:  Configuration block to read parameter from.
 * param: Name of parameter to load.
 */
slibreal_t InitAvalancheDistribution_param(ConfigBlock *conf, const string &param) {
    Setting *set;
    if (!conf->HasSetting(param))
        throw SOFTException("Distribution function '%s': Required parameter '%s' not specified.", conf->GetName().c_str(), param.c_str());

    set = conf->GetSetting(param);
    if (!set->IsScalar())
        throw SOFTException("Distribution function '%s': %s: Invalid value assigned to parameter. Expected real scalar.", conf->GetName().c_str(), param.c_str());

    return set->GetScalar();
}

/**
 * Initialize distribution function from CODE output.
 */
RadialDistributionFunction *InitCODEDistribution(MagneticField2D *magfield, ConfigBlock *conf, ConfigBlock *root) {
    RadialProfile *radprof;
    Setting *set;
    string name;
    int interptype = CODEDistributionFunction::INTERPOLATION_CSPLINE;
    int timestep = -1;
    
    // Name
    if (!conf->HasSetting("name"))
        throw SOFTException("Distribution function '%s': Name of numerical distribution function file not specified.", conf->GetName().c_str());
    else {
        set = conf->GetSetting("name");
        if (set->GetNumberOfValues() != 1)
            throw SOFTException("Distribution function '%s': name: Invalid value assigned to parameter. Expected string.", conf->GetName().c_str());

        name = set->GetString();
    }

    // Type of interpolation
    if (conf->HasSetting("interptype")) {
        set = conf->GetSetting("interptype");
        if (set->GetNumberOfValues() != 1)
            throw SOFTException("Distribution function '%s': interptype: Invalid value assigned to parameter. Expected string.", conf->GetName().c_str());

        if (set->GetString() == "linear")
            interptype = CODEDistributionFunction::INTERPOLATION_LINEAR;
        else if (set->GetString() == "polynomial")
            interptype = CODEDistributionFunction::INTERPOLATION_POLYNOMIAL;
        else if (set->GetString() == "cspline")
            interptype = CODEDistributionFunction::INTERPOLATION_CSPLINE;
        else if (set->GetString() == "cspline_periodic")
            interptype = CODEDistributionFunction::INTERPOLATION_CSPLINE_PERIODIC;
        else if (set->GetString() == "akima")
            interptype = CODEDistributionFunction::INTERPOLATION_AKIMA;
        else if (set->GetString() == "akima_periodic")
            interptype = CODEDistributionFunction::INTERPOLATION_AKIMA_PERIODIC;
        else if (set->GetString() == "steffen")
            interptype = CODEDistributionFunction::INTERPOLATION_STEFFEN;
        else
            throw SOFTException(
                "Distribution function '%s': interptype: Unrecognized interpolation type requested: %s.",
                conf->GetName().c_str(), set->GetString().c_str()
            );
    }

    // Timestep
    if (conf->HasSetting("time")) {
        set = conf->GetSetting("time");
        if (set->IsInteger32())
            timestep = set->GetInteger32();
        else if (set->GetNumberOfValues() == 1 && 
            (set->GetString() == "end" || set->GetString() == "last"))
            timestep = -1;
        else
            throw SOFTException("Distribution function '%s': time: Invalud value assigned to parameter. Expected 32-bit integer.", conf->GetName().c_str());
    }

    // Load radial profile (if specified)
    if (conf->HasSetting("radprof")) {
        set = conf->GetSetting("radprof");
        radprof = InitRadialProfile(magfield, set, root, conf->GetName());
    } else
        radprof = new UniformRadialProfile();

    // Combine...
    return new RadialDistributionFunction(
        radprof,
        new CODEDistributionFunction(
            name, timestep, interptype
        )
    );
}

/**
 * Initialize distribution function from LUKE output.
 */
DistributionFunction *InitLUKEDistribution(MagneticField2D *magfield, ConfigBlock *conf) {
    Setting *set;
    string name;
    bool logarithmize = false;
    int interptype = NumericMomentumSpaceDistributionFunction::INTERPOLATION_CUBIC;
    
    // Name
    if (!conf->HasSetting("name"))
        throw SOFTException("Distribution function '%s': Name of numerical distribution function file not specified.", conf->GetName().c_str());
    else {
        set = conf->GetSetting("name");
        if (set->GetNumberOfValues() != 1)
            throw SOFTException("Distribution function '%s': name: Invalid value assigned to parameter. Expected string.", conf->GetName().c_str());

        name = set->GetString();
    }

    // Interpolate logarithmically
    if (conf->HasSetting("logarithmize")) {
        set = conf->GetSetting("logarithmize");
        if (!set->IsBool())
            throw SOFTException("Distribution function '%s': logarithmize: Invalid value assigned to parameter. Expected boolean.", conf->GetName().c_str());

        logarithmize = set->GetBool();
    }

    // Type of interpolation
    if (conf->HasSetting("interptype")) {
        set = conf->GetSetting("interptype");
        if (set->GetNumberOfValues() != 1)
            throw SOFTException("Distribution function '%s': interptype: Invalid value assigned to parameter. Expected string.", conf->GetName().c_str());

        if (set->GetString() == "cubic")
            interptype = NumericMomentumSpaceDistributionFunction::INTERPOLATION_CUBIC;
        else if (set->GetString() == "linear")
            interptype = NumericMomentumSpaceDistributionFunction::INTERPOLATION_LINEAR;
        else
            throw SOFTException(
                "Distribution function '%s': interptype: Unrecognized interpolation type requested: %s.",
                conf->GetName().c_str(), set->GetString().c_str()
            );
    }

    return new LUKEDistributionFunction(name, magfield, logarithmize, interptype);
}

/**
 * Initialize distribution function from NORSE output.
 */
DistributionFunction *InitNORSEDistribution(MagneticField2D*, ConfigBlock *conf, ConfigBlock*) {
    throw SOFTException("Distribution function '%s': type: Support for NORSE distribution functions has not been implemented yet.", conf->GetName().c_str());
}

/**
 * Initialize distribution function from a 3D numerical distribution
 * (aka SOFT distribution).
 */
SOFTDistributionFunction *InitNumericalDistribution(MagneticField2D *magfield, ConfigBlock *conf) {
    Setting *set;
    string name;
    bool logarithmize = false;
    int interptype = NumericMomentumSpaceDistributionFunction::INTERPOLATION_CUBIC;
    
    // Name
    if (!conf->HasSetting("name"))
        throw SOFTException("Distribution function '%s': Name of numerical distribution function file not specified.", conf->GetName().c_str());
    else {
        set = conf->GetSetting("name");
        if (set->GetNumberOfValues() != 1)
            throw SOFTException("Distribution function '%s': name: Invalid value assigned to parameter. Expected string.", conf->GetName().c_str());

        name = set->GetString();
    }

    // Interpolate logarithmically
    if (conf->HasSetting("logarithmize")) {
        set = conf->GetSetting("logarithmize");
        if (!set->IsBool())
            throw SOFTException("Distribution function '%s': logarithmize: Invalid value assigned to parameter. Expected boolean.", conf->GetName().c_str());

        logarithmize = set->GetBool();
    }

    // Type of interpolation
    if (conf->HasSetting("interptype")) {
        set = conf->GetSetting("interptype");
        if (set->GetNumberOfValues() != 1)
            throw SOFTException("Distribution function '%s': interptype: Invalid value assigned to parameter. Expected string.", conf->GetName().c_str());

        if (set->GetString() == "cubic")
            interptype = NumericMomentumSpaceDistributionFunction::INTERPOLATION_CUBIC;
        else if (set->GetString() == "linear")
            interptype = NumericMomentumSpaceDistributionFunction::INTERPOLATION_LINEAR;
        else
            throw SOFTException(
                "Distribution function '%s': interptype: Unrecognized interpolation type requested: %s.",
                conf->GetName().c_str(), set->GetString().c_str()
            );
    }

    return new SOFTDistributionFunction(name, magfield, logarithmize, interptype);
}

/**
 * Initialize an exponential pitch distribution function.
 *
 * magfield: Magnetic field that the distribution lives in.
 * conf:     Configuration of module.
 */
RadialDistributionFunction *InitExponentialPitchDistribution(MagneticField2D *magfield, ConfigBlock *conf, ConfigBlock *root) {
    slibreal_t C;

    if (!conf->HasSetting("C")) {
        throw SOFTException(
            "Distribution function '%s': Pitch distribution exponent 'C' not specified.",
            conf->GetName().c_str()
        );
    } else {
        Setting *s = conf->GetSetting("C");
        if (!s->IsScalar())
            throw SOFTException(
                "Distribution function '%s': The given pitch distribution exponent 'C' is not a scalar.",
                conf->GetName().c_str()
            );

        C = s->GetScalar();
    }

    RadialProfile *radprof;
    if (conf->HasSetting("radprof")) {
        Setting *set = conf->GetSetting("radprof");
        radprof = InitRadialProfile(magfield, set, root, conf->GetName());
    } else
        radprof = new UniformRadialProfile();

    // Combine...
    return new RadialDistributionFunction(
        radprof, new ExponentialPitch(C)
    );
}

RadialDistributionFunction *InitUnitDistributionFunction(MagneticField2D *magfield, ConfigBlock *conf, ConfigBlock *root) {
    RadialProfile *radprof;

    // Load radial profile (if specified)
    if (conf->HasSetting("radprof")) {
        Setting *set = conf->GetSetting("radprof");
        radprof = InitRadialProfile(magfield, set, root, conf->GetName());
    } else
        radprof = new UniformRadialProfile();

    // Combine...
    return new RadialDistributionFunction(
        radprof, new UnitDistributionFunction()
    );
}

