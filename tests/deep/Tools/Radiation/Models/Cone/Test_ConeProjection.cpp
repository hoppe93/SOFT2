/**
 * Test various cone models. Tests included:
 *
 * [X] Compare output of 'Original' model to the model
 *     implemented in SOFTv1.
 */

#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/RadiationParticle.h"
#include "Tools/Radiation/Models/Cone.h"
#include "Tools/Radiation/Models/Cone/Projection/ConeProjection.h"
#include "Tools/Radiation/Models/Cone/Projection/Original.h"
#include "Tools/Radiation/Models/Cone/Projection/Reverse.h"
#include "SOFTException.h"
#include "Test_ConeProjection.h"

// Include lookup table
#include "Test_ConeProjection_tbls.cpp"

using namespace __Radiation;

/**
 * Compare the output of the given model to the
 * tabulated output of the SOFTv1 cone model.
 *
 * cp:       Cone projection model to test.
 * tabindex: Index of lookup table to compare to.
 * tol:      Allowed tolerance.
 * name:     Name of the projection model tested.
 */
bool Test_ConeProjection::CompareModelToTable(ConeProjection *cp, unsigned int tabindex, slibreal_t tol, string &name) {
    unsigned int i;
    slibreal_t fraction, csfraction, Delta;

    for (i = 0; i < conespec_N[tabindex]; i++) {
        // Generate RadiationParticle
        RadiationParticle *rp = GetRadiationParticle(conespec[tabindex], i, cp->GetDetector());
        csfraction = conespec[tabindex][i].fraction;

        // Evaluate SOFTv2 cone model
        fraction = cp->ComputeOverlappingFraction(rp);
        Delta = fabs((fraction - csfraction) / csfraction);

        delete rp;

        if (Delta > tol) {
            this->PrintError("Cone model '%s' failed at index [%u.%u] Delta = %e.", name.c_str(), tabindex, i, Delta);
            return false;
        }// else printf("Delta = %e\n", Delta);
    }

    return true;
}

/**
 * Generate a 'RadiationParticle' object from the
 * entry in 'conespec' with index 'index'.
 *
 * index: Index of the entry in the 'conespec' lookup
 *        table to generate the RadiationParticle from.
 *
 * RETURNS a new RadiationParticle object which must be
 * deallocated after use.
 */
RadiationParticle *Test_ConeProjection::GetRadiationParticle(const struct conemodel_spec *tbl, unsigned int index, Detector *det) {
    const struct conemodel_spec *cs = tbl+index;
    slibreal_t Jdtdrho = 1.0,
               ppar = cs->ppar,
               pperp = cs->pperp,
               p2 = ppar*ppar + pperp*pperp,
               gamma = sqrt(1 + p2),
               B = 1.0,
               m = ELECTRON_MASS,
               q = -ELEMENTARY_CHARGE,
               Bvec[3] = {1.0,0.0,0.0};

    Vector<3> p(cs->vhat);
    // The magnitude of p doesn't really matter,
    // so we set it approximately to ppar (~ GC momentum)
    p = ppar * p;

    return new RadiationParticle(
        cs->x, p, Jdtdrho, ppar, pperp,
        gamma, p2, det->GetPosition(), B,
        Bvec, Bvec, m, q, 0, 0, 0
    );
}

/**
 * [TEST] compare cone models to the cone
 * model in SOFTv1 using table lookup.
 */
bool Test_ConeProjection::RunSOFTv1Comparison() {
    bool success = true;

    // Original model
    success &= CompareToAllTables<ConeProjectionOriginal>("Original", 2e-6);

    // Reverse model
    success &= CompareToAllTables<ConeProjectionReverse>("Reverse", 1e-3);
    
    return success;
}

template<class T>
bool Test_ConeProjection::CompareToAllTables(string name, slibreal_t tol) {
    bool success = true, psuccess = true;
    unsigned int i;
    slibreal_t tilt = 0;

    for (i = 0; i < Test_ConeProjection::nconespec; i++) {
        Vector<3> dir(detdir[i]), pos(detpos[i]);
        Detector *det = new Detector(
            aperture[i], tilt, visang[i],
            dir, pos, 0
        );

        T *cpo = new T(det);
        psuccess = CompareModelToTable(cpo, i, tol, name);

        if (psuccess)
            this->PrintOK("Model '%s' compared successfully to SOFTv1 table #%u.", name.c_str(), (i+1));

        success &= psuccess;

        delete cpo;
        delete det;
    }

    return success;
}

/**
 * Run all tests of this module.
 */
bool Test_ConeProjection::Run(bool) {
    bool success = true;

    // Compare 'Original' model to SOFTv1 model
    success &= RunSOFTv1Comparison();

    return success;
}

