#ifndef _RADIATION_PARTICLE_H
#define _RADIATION_PARTICLE_H

#include "Orbit/Orbit.h"

namespace __Radiation {
    class RadiationParticle {
        private:
            Vector<3> x, p, phat, rcp, rcphat, Bvec, bHat, gradB, curlB, curlBhat;
            slibreal_t **jacobianB=nullptr;
            slibreal_t B, R, diffel;
            slibreal_t ppar, pperp, gamma, p2;
            slibreal_t dphi, Rdphi, Jdtdrho, drho;
            slibreal_t dp1, dp2;

            slibreal_t rcp2, rcplen;
            slibreal_t m, q;            // Particle mass and charge
            slibreal_t f;               // Distribution function value
            
            unsigned int ir, ip1, ip2;
        public:
            RadiationParticle(Orbit*, unsigned int, Vector<3>&);
            RadiationParticle(
                const Vector<3>&, const Vector<3>&, slibreal_t, slibreal_t,
                slibreal_t, slibreal_t, slibreal_t, const Vector<3>&,
                slibreal_t, const Vector<3>&, const Vector<3>&,
                slibreal_t, slibreal_t, unsigned int, unsigned int, unsigned int,
                slibreal_t *gradB=nullptr, slibreal_t *curlB=nullptr, slibreal_t **jacobianB=nullptr
            );
            virtual ~RadiationParticle();

            slibreal_t GetF()                   const { return f; }       // Distribution function value
            Vector<3>& GetPosition()                  { return x; }      // Particle position
            Vector<3>& GetRCP()                       { return rcp; }         // Detector-particle vector
            Vector<3>& GetRCPHat()                    { return rcphat; }   // Normalized detector-particle vector
            Vector<3>& GetP()                         { return p; }
            Vector<3>& GetPHat()                      { return phat; }
            slibreal_t GetRCP2()                const { return rcp2; }       // Length-squared of detector-particle vector
            slibreal_t GetRCPLength()           const { return rcplen; }// Length of detector-particle vector
            slibreal_t GetDphi()                const { return dphi; }       // Toroidal step (angular difference)
            slibreal_t GetRDphi()               const { return Rdphi; }     // Radial (cylindrical) particle distance times dphi
            slibreal_t GetJdtdrho()             const { return Jdtdrho; } // Spatial Jacobian (times drho and dtau)
            slibreal_t GetDrho()                const { return drho; }       // Radial step
            slibreal_t GetDifferentialElement() const { return diffel; }  // "Total" real+momentum space differential element

            slibreal_t GetB()                   const { return B; }
            Vector<3>& GetBvec()                      { return Bvec; }
            Vector<3>& GetBhat()                      { return bHat; }
            slibreal_t GetGamma()               const { return gamma; }
            slibreal_t GetPpar()                const { return ppar; }
            slibreal_t GetPperp()               const { return pperp; }
            slibreal_t GetP2()                  const { return p2; }

            slibreal_t GetDP1()                 const { return dp1; }
            slibreal_t GetDP2()                 const { return dp2; }

            slibreal_t GetCharge()              const { return q; }
            slibreal_t GetMass()                const { return m; }

            Vector<3>& GetGradB()                     { return gradB; }
            Vector<3>& GetCurlB()                     { return curlB; }
            Vector<3>& GetCurlBHat()                  { return curlBhat; }
            slibreal_t **GetBJacobian()         const { return jacobianB; }

            unsigned int GetIndexR()            const { return ir; }
            unsigned int GetIndexP1()           const { return ip1; }
            unsigned int GetIndexP2()           const { return ip2; }

            void SetDistributionValue(const slibreal_t f) { this->f = f; }
            void SetDifferentialElements(slibreal_t, slibreal_t, slibreal_t, slibreal_t);
            void UpdateXY(slibreal_t, slibreal_t, slibreal_t, slibreal_t, Vector<3>&);
    };
}

#endif/*_RADIATION_PARTICLE_H*/
