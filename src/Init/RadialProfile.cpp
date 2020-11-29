/**
 * Initialization of the radial profile associated with
 * the distribution function.
 */

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/DistributionFunction/BesselRadialProfile.h>
#include <softlib/DistributionFunction/LinearRadialProfile.h>
#include <softlib/DistributionFunction/PowerRadialProfile.h>
#include <softlib/DistributionFunction/GaussianRadialProfile.h>
#include <softlib/DistributionFunction/UniformRadialProfile.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include "Init/InitConfig.h"
#include "Init/InitDistribution.h"
#include "SOFTException.h"

using namespace std;

/**
 * Initialize a radial profile that can be coupled to a
 * momentum space distribution function.
 */
RadialProfile *InitRadialProfile(MagneticField2D *magfield, Setting *rpset, ConfigBlock *root, const string &dfname) {
    string type;

    if (rpset->GetNumberOfValues() != 1)
        throw SOFTException("Distribution function '%s': radprof: Invalid value assigned to parameter. Expected string.", dfname.c_str());

    string name = rpset->GetString();
    if (!root->HasSubBlock(CONFBLOCK_RADIALPROFILE, name)) {
        if (name == "uniform")
            return new UniformRadialProfile();
        else if (name == "linear")
            return InitLinearRadialProfile(magfield);
        else
            throw SOFTException("Distribution function '%s': No radial profile named '%s' defined in the configuration.", dfname.c_str(), name.c_str());
    }

    ConfigBlock *conf = root->GetConfigBlock(CONFBLOCK_RADIALPROFILE, name);

    type = conf->GetSecondaryType();
    if (type == "bessel")
        return InitBesselRadialProfile(magfield, conf);
    else if (type == "linear")
        return InitLinearRadialProfile(magfield, conf);
    else if (type == "power")
        return InitPowerRadialProfile(magfield, conf);
    else if (type == "gaussian")
        return InitGaussianRadialProfile(magfield, conf);
    else if (type == "uniform")
        return new UniformRadialProfile();
    else
        throw SOFTException("Radial profile '%s': type: Unrecognized value assigned to parameter.", name.c_str());
}

/**
 * Initialize a Bessel J0 radial profile using the
 * given configuration block.
 */
BesselRadialProfile *InitBesselRadialProfile(MagneticField2D *magfield, ConfigBlock *conf) {
    slibreal_t rmin, rmax;

    rmin = magfield->GetMagneticAxisR();
    rmax = magfield->GetMaxRadius();

    if (conf != nullptr) {
        InitRadialProfile_get_radial_limits(
            magfield, conf, &rmin, &rmax
        );
    }

    return new BesselRadialProfile(rmin, rmax);
}

/**
 * Initialize a linear radial profile using the
 * given configuration block.
 */
LinearRadialProfile *InitLinearRadialProfile(MagneticField2D *magfield, ConfigBlock *conf) {
    slibreal_t rmin, rmax;

    rmin = magfield->GetMagneticAxisR();
    rmax = magfield->GetMaxRadius();

    if (conf != nullptr) {
        InitRadialProfile_get_radial_limits(
            magfield, conf, &rmin, &rmax
        );
    }

    return new LinearRadialProfile(rmin, rmax);
}

/**
 * Initialize a power-function radial profile using the
 * given configuration block.
 */
PowerRadialProfile *InitPowerRadialProfile(MagneticField2D *magfield, ConfigBlock *conf) {
    Setting *set;
    slibreal_t rmin, rmax, b;

    InitRadialProfile_get_radial_limits(
        magfield, conf, &rmin, &rmax
    );

    if (!conf->HasSetting("exponent"))
        throw SOFTException("Radial profile '%s': exponent: No value assigned to mandatory parameter.", conf->GetName().c_str());

    set = conf->GetSetting("exponent");
    if (!set->IsScalar())
        throw SOFTException("Radial profile '%s': exponent: Invalid value assigned to parameter.", conf->GetName().c_str());

    b = set->GetScalar();

    return new PowerRadialProfile(rmin, rmax, b);
}
/**
 * Initialize a gaussian-function radial profile using the
 * given configuration block.
 */
GaussianRadialProfile *InitGaussianRadialProfile(MagneticField2D *magfield, ConfigBlock *conf) {
    Setting *set_a,*set_b,*set_c;
    slibreal_t rmin, rmax, a, b, c;

    InitRadialProfile_get_radial_limits(
        magfield, conf, &rmin, &rmax
    );

    if (!conf->HasSetting("gau_a"))
        throw SOFTException("Radial profile '%s': first-a param: No value assigned to mandatory parameter.", conf->GetName().c_str());
    if (!conf->HasSetting("gau_b"))
        throw SOFTException("Radial profile '%s': second-b param: No value assigned to mandatory parameter.", conf->GetName().c_str());
    if (!conf->HasSetting("gau_c"))
        throw SOFTException("Radial profile '%s': third-c param: No value assigned to mandatory parameter.", conf->GetName().c_str());

    set_a = conf->GetSetting("gau_a");
    if (!set_a->IsScalar())
        throw SOFTException("Radial profile '%s': exponent: Invalid value assigned to parameter.", conf->GetName().c_str());

    set_b = conf->GetSetting("gau_b");
    if (!set_b->IsScalar())
        throw SOFTException("Radial profile '%s': exponent: Invalid value assigned to parameter.", conf->GetName().c_str());

    set_c = conf->GetSetting("gau_c");
    if (!set_c->IsScalar())
        throw SOFTException("Radial profile '%s': exponent: Invalid value assigned to parameter.", conf->GetName().c_str());

    a = set_a->GetScalar();
    b = set_b->GetScalar();
    c = set_c->GetScalar();

    return new GaussianRadialProfile(rmin, rmax, a, b, c);
}
/**
 * Read the limits of a radial coordinate from the
 * given configuration object. This function allows
 * the specification of either the lower (*min) or
 * upper (*max) radial coordinate, as either the major
 * radius (rho*), minor radius (r*) or normalized
 * minor radius (a*).
 */
void InitRadialProfile_get_radial_limits(
    MagneticField2D *magfield, ConfigBlock *conf,
    slibreal_t *r0, slibreal_t *r1
) {
    slibreal_t rmin, rmax;
    *r0 = rmin = magfield->GetMagneticAxisR();
    *r1 = rmax = magfield->GetMaxRadius();

    if (conf == nullptr)
        return;

    // Inner radius
    if (conf->HasSetting("rhomin"))
        *r0 = InitRadialProfile_get_radial_limits_inner(conf->GetSetting("rhomin"), rmin, rmax, conf->GetName());
    else if (conf->HasSetting("rmin"))
        *r0 = rmin + InitRadialProfile_get_radial_limits_inner(conf->GetSetting("rmin"), 0.0, rmax-rmin, conf->GetName());
    else if (conf->HasSetting("amin"))
        *r0 = rmin + (rmax-rmin)*InitRadialProfile_get_radial_limits_inner(conf->GetSetting("amin"), 0.0, 1.0, conf->GetName());
    else
        *r0 = rmin;

    // Outer radius
    if (conf->HasSetting("rhomax"))
        *r1 = InitRadialProfile_get_radial_limits_inner(conf->GetSetting("rhomax"), rmin, rmax, conf->GetName());
    else if (conf->HasSetting("rmax"))
        *r1 = rmin + InitRadialProfile_get_radial_limits_inner(conf->GetSetting("rmax"), 0.0, rmax-rmin, conf->GetName());
    else if (conf->HasSetting("amax"))
        *r1 = rmin + (rmax-rmin)*InitRadialProfile_get_radial_limits_inner(conf->GetSetting("amax"), 0.0, 1.0, conf->GetName());
    else
        *r1 = rmax;

    if (*r0 >= rmax)
        throw SOFTException("Radial profile '%s': Inner radius is larger than outer radius.", conf->GetName().c_str());
}

/**
 * Loads the value of the given setting, applying the
 * typical checks for radial coordinates.
 *
 * set:    Settings object specifying the setting.
 * lowlim: Lowest permitted value.
 * uplim:  Highest permitted value.
 */
slibreal_t InitRadialProfile_get_radial_limits_inner(
    Setting *set, const slibreal_t lowlim, const slibreal_t, const string &rpname
) {
    slibreal_t tmp;
    if (!set->IsScalar())
        throw SOFTException(
            "Radial profile '%s': %s: Invalid value assigned to parameter. Expected scalar.",
            rpname.c_str(), set->GetName().c_str()
        );

    tmp = set->GetScalar();
    if (tmp < lowlim) {
        throw SOFTException(
            "Radial profile '%s': %s: Invalid value assigned to parameter. Lower point must be greater than %f.",
            rpname.c_str(), set->GetName().c_str(), lowlim
        );
    /*} else if (tmp > uplim) {
        SOFT::PrintWarning(...);*/
    }

    return tmp;
}
