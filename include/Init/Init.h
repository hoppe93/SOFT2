#ifndef _INIT_SOFT_H
#define _INIT_SOFT_H

#include <string>
#include <softlib/Configuration.h>

#include "SOFT.h"

SOFT *InitSOFT(Configuration*);
std::string InitSOFTDefaults();
slibreal_t init_get_scalar(ConfigBlock*, const std::string&, const std::string&);
uint32_t init_get_uint32(ConfigBlock*, const std::string&, const std::string&);
std::string& init_get_string(ConfigBlock*, const std::string&, const std::string&);

#include "Init/InitConfig.h"
#include "Init/InitDistribution.h"

#endif/*_INIT_SOFT_H*/
