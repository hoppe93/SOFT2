#ifndef _PROJECTION_ORIGINAL_H
#define _PROJECTION_ORIGINAL_H

#include <softlib/Vector.h>
#include "Tools/Radiation/Models/Cone/Projection/ConeProjection.h"

namespace __Radiation {
    struct intersection_points {
        slibreal_t t[8], x[8], y[8];
        unsigned int nt;
    };

    class ConeProjectionOriginal : public ConeProjection {
        private:
            Vector<3> e1r, e2r, Pne1r, Pne2r;
            slibreal_t Pn[3][3];
            slibreal_t absPne1r, absPne2r;
        public:
            ConeProjectionOriginal(Detector *det) : ConeProjection(det) {}

            virtual slibreal_t ComputeOverlappingFraction(RadiationParticle*);

            void EllipseGetIntersections(
                slibreal_t, slibreal_t, slibreal_t,
                slibreal_t, slibreal_t,
                slibreal_t, slibreal_t,
                struct intersection_points*
            );
            void EllipseGetIntersection(
                slibreal_t, slibreal_t, slibreal_t,
                slibreal_t, slibreal_t, slibreal_t,
                slibreal_t, slibreal_t[8], slibreal_t[8],
                struct intersection_points*
            );
            void HyperbolaGetIntersections(
                slibreal_t, slibreal_t, slibreal_t,
                slibreal_t, slibreal_t,
                slibreal_t, slibreal_t,
                struct intersection_points*
            );
            void HyperbolaGetIntersection(
                slibreal_t, slibreal_t, slibreal_t,
                slibreal_t, slibreal_t, slibreal_t,
                slibreal_t, slibreal_t[8], slibreal_t[8],
                struct intersection_points*
            );
            void InitReverseProjection(Vector<3>&, slibreal_t, slibreal_t);
            void SortIntersectionPoints(struct intersection_points&);
            slibreal_t ToPointOnUnitCone(slibreal_t, slibreal_t);
    };
};

#endif/*_PROJECTION_ORIGINAL_H*/
