#ifndef _PROJECTION_REVERSE_H
#define _PROJECTION_REVERSE_H

#include <softlib/Vector.h>
#include "Tools/Radiation/Models/Cone/Projection/ConeProjection.h"

namespace __Radiation {
    class ConeProjectionReverse : public ConeProjection {
        private:
        public:
            ConeProjectionReverse(Detector *det) : ConeProjection(det) {}

            virtual slibreal_t ComputeOverlappingFraction(RadiationParticle*);
            void FindCircleIntersectionPoints(
                Vector<3>&, Vector<3>&, Vector<3>&,
                Vector<3>&, Vector<3>&, slibreal_t,
                unsigned int&, slibreal_t[8], unsigned int&
            );
            slibreal_t GetParameterT(const slibreal_t, const slibreal_t, const slibreal_t);
            void SetupBasisVectors(Vector<3>&, Vector<3>&, Vector<3>&);
    };
};

#endif/*_PROJECTION_REVERSE_H*/
