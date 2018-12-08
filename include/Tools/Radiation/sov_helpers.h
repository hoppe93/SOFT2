#ifndef _RADIATION_SOV_HELPERS_H
#define _RADIATION_SOV_HELPERS_H

#include <functional>
#include "config.h"

namespace __Radiation {
    void bracket_min(std::function<slibreal_t(const slibreal_t)>, const slibreal_t, slibreal_t*, slibreal_t*);
    bool bracket_root(std::function<slibreal_t(const slibreal_t)>, const slibreal_t, const slibreal_t, const unsigned int, slibreal_t*, slibreal_t*);

    slibreal_t find_min(std::function<slibreal_t(const slibreal_t)>, const slibreal_t, const slibreal_t);
    slibreal_t find_root(std::function<slibreal_t(const slibreal_t)>, slibreal_t, slibreal_t, const slibreal_t);

    slibreal_t mod2pi(const slibreal_t);
    void shift3(slibreal_t&, slibreal_t&, slibreal_t&, const slibreal_t);
}

#endif/*_RADIATION_SOV_HELPERS_H*/
