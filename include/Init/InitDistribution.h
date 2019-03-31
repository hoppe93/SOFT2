#ifndef _INIT_DISTRIBUTION_H
#define _INIT_DISTRIBUTION_H

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/DistributionFunction/DistributionFunction.h>
#include <softlib/DistributionFunction/LinearRadialProfile.h>
#include <softlib/DistributionFunction/PowerRadialProfile.h>
#include <softlib/DistributionFunction/RadialProfile.h>
#include <softlib/DistributionFunction/RadialDistributionFunction.h>
#include <softlib/DistributionFunction/SOFTDistributionFunction.h>
#include <softlib/MagneticField/MagneticField2D.h>

DistributionFunction *InitDistributionFunction(MagneticField2D*, ConfigBlock*, ConfigBlock*);
RadialDistributionFunction *InitAvalancheDistribution(MagneticField2D*, ConfigBlock*, ConfigBlock*);
slibreal_t InitAvalancheDistribution_param(ConfigBlock*, const std::string&);
RadialDistributionFunction *InitCODEDistribution(MagneticField2D*, ConfigBlock*, ConfigBlock*);
DistributionFunction *InitLUKEDistribution(MagneticField2D*, ConfigBlock*);
DistributionFunction *InitNORSEDistribution(MagneticField2D*, ConfigBlock*, ConfigBlock*);
SOFTDistributionFunction *InitNumericalDistribution(MagneticField2D*, ConfigBlock*);
RadialDistributionFunction *InitUnitDistributionFunction(MagneticField2D*, ConfigBlock*, ConfigBlock*);

// Radial profiles
RadialProfile *InitRadialProfile(MagneticField2D*, Setting*, ConfigBlock*, const std::string&);
LinearRadialProfile *InitLinearRadialProfile(MagneticField2D*, ConfigBlock *conf=nullptr);
PowerRadialProfile *InitPowerRadialProfile(MagneticField2D*, ConfigBlock*);
void InitRadialProfile_get_radial_limits(MagneticField2D*, ConfigBlock*, slibreal_t*, slibreal_t*);
slibreal_t InitRadialProfile_get_radial_limits_inner(Setting*, const slibreal_t, const slibreal_t, const std::string&);

#endif/*_INIT_DISTRIBUTION_H*/
