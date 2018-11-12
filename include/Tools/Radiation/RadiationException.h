#ifndef _RADIATION_EXCEPTION_H
#define _RADIATION_EXCEPTION_H

#include <string>
#include <utility>
#include "Tools/Tool.h"

namespace __Radiation {
    class RadiationException : public ToolException {
        public:
            template<typename ... Args>
            RadiationException(const std::string &msg, Args&& ... args)
                : ToolException(msg, std::forward<Args>(args) ...) {
                AddModule("Radiation");
            }
    };
}

#endif/*_RADIATION_EXCEPTION_H*/
