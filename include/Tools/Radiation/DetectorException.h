#ifndef _DETECTOR_EXCEPTION_H
#define _DETECTOR_EXCEPTION_H

#include "RadiationException.h"

namespace __Radiation {
    class DetectorException : public RadiationException {
        public:
            template<typename ... Args>
            DetectorException(const std::string &msg, Args&& ... args)
                : RadiationException(msg, std::forward<Args>(args) ...) {
                AddModule("Detector");
            }
    };
}

#endif/*_DETECTOR_EXCEPTION_H*/
