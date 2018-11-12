#ifndef _ANGULAR_DISTRIBUTION_EXCEPTION_H
#define _ANGULAR_DISTRIBUTION_EXCEPTION_H

#include <string>
#include "Tools/Radiation/RadiationException.h"

namespace __Radiation {
    class AngularDistributionException : public RadiationException {
        public:
            template<typename ... Args>
            AngularDistributionException(const std::string &msg, Args&& ... args)
                : RadiationException(msg, std::forward<Args>(args) ...) {
                AddModule("Angular distribution");
            }
    };
}

#endif/*_ANGULAR_DISTRIBUTION_EXCEPTION_H*/
