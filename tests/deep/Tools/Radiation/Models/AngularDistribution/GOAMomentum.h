#ifndef _GOA_MOMENTUM_H
#define _GOA_MOMENTUM_H

#include <softlib/config.h>
#include <softlib/Vector.h>
#include <softlib/MagneticField/MagneticField2D.h>

class GOAMomentum {
    private:
        Vector<3> B, bhat;
        Vector<3> gradB, curlB;
        Vector<3> curlBhat, kappa;
        Vector<3> one, two;
        slibreal_t J[3][3];
        slibreal_t Babs;
        slibreal_t charge;

        MagneticField2D *magfield;
    public:
        GOAMomentum(Vector<3>&, MagneticField2D*, slibreal_t);

        Vector<3>& Acceleration(Vector<3>&, slibreal_t, slibreal_t, slibreal_t, Vector<3>&);
        void GetContractions(Vector<3>&, Vector<3>&, slibreal_t*, slibreal_t*);
        Vector<3>& GuidingCenterMomentum(slibreal_t, slibreal_t, Vector<3>&);
        Vector<3>& GyroVector(slibreal_t, slibreal_t, slibreal_t, Vector<3>&); 
        void InitLocalCoordinateSystem();
        Vector<3>& Momentum1(slibreal_t, slibreal_t, slibreal_t, Vector<3>&);
        Vector<3>& ParticleMomentum(slibreal_t, slibreal_t, slibreal_t, Vector<3>&, Vector<3>&);
        void SetGCPosition(Vector<3>&, MagneticField2D*);

        slibreal_t GetBabs() { return Babs; }
        Vector<3>& GetBhat() { return bhat; }
        Vector<3>& GetOne()  { return one; }
        Vector<3>& GetTwo()  { return two; }
        void DumpMF();
};
#endif/*_GOA_MOMENTUM_H*/
