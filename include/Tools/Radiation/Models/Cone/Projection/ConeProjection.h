#ifndef _CONE_PROJECTION_H
#define _CONE_PROJECTION_H

#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"

namespace __Radiation {
    class ConeProjection {
        protected:
            Detector *detector;
        public:
            ConeProjection(Detector *det) { this->detector = det; }
            virtual ~ConeProjection() {};
            virtual slibreal_t ComputeOverlappingFraction(RadiationParticle*) = 0;

            Detector *GetDetector() { return detector; }
    };
}
#endif/*_CONE_PROJECTION_H*/
