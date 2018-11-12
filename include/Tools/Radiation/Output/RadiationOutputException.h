#ifndef _RADIATION_OUTPUT_EXCEPTION_H
#define _RADIATION_OUTPUT_EXCEPTION_H

#include "Tools/Radiation/Radiation.h"
#include "Tools/Radiation/RadiationException.h"

namespace __Radiation {
    class RadiationOutputException : public RadiationException {
        public:
            template<typename ... Args>
            RadiationOutputException(const std::string &msg, Args&& ... args)
                : RadiationException(msg, std::forward<Args>(args) ...) {
                AddModule("Output");
            }
    };
}

#endif/*_RADIATION_OUTPUT_EXCEPTION_H*/
