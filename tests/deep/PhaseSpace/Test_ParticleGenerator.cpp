/**
 * Unit test for the 'ParticleGenerator' class.
 * The following tests are carried out:
 *
 * - Check that the phase-space grid is properly generated.
 * - TODO Check that the effective magnetic axis is properly calculated.
 */

#include <string>

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include "Orbit/GuidingCenterEquation.h"
#include "PhaseSpace/ParticleGenerator.h"

#include "Test_ParticleGenerator.h"

using namespace std;

const slibreal_t Test_ParticleGenerator::effax_ppar0    = 5.0;
const slibreal_t Test_ParticleGenerator::effax_ppar1    = 100.0;
const slibreal_t Test_ParticleGenerator::effax_pperp0   = 5.0;
const slibreal_t Test_ParticleGenerator::effax_pperp1   = 100.0;
const unsigned int Test_ParticleGenerator::effax_nppar  = 20;
const unsigned int Test_ParticleGenerator::effax_npperp = 20;

/**
 * Calculates the guiding-center drift shift using the
 * 'GuidingCenterEquation' implementation of the
 * guiding-center equations of motion.
 */
slibreal_t Test_ParticleGenerator::GetDriftShift(
    MagneticField2D *mf, GuidingCenterEquation *eqn,
    const slibreal_t ppar, const slibreal_t pperp,
    const slibreal_t tol
) {
    // Lambda that evaluates Xpol-dot (poloidal
    // guiding-center speed
    auto gxpol = [&ppar,&pperp,&mf,&eqn](slibreal_t x) {
        Vector<6> zval, dzdt;
        Vector<3> Xdot, B, bhat;
        Particle p;
        slibreal_t
            X0 = x,
            Z0 = mf->GetMagneticAxisZ(),
            c = LIGHTSPEED,
            m = ELECTRON_MASS,
            q = -ELEMENTARY_CHARGE,
            Babs;

        B = mf->Eval(X0, 0.0, Z0);
        Babs = B.Norm();
        bhat = B / Babs;

        zval[GuidingCenterEquation::COORD_Y]    = 0.0;
        zval[GuidingCenterEquation::COORD_Z]    = Z0;
        zval[GuidingCenterEquation::COORD_PPAR] = ppar;
        zval[GuidingCenterEquation::COORD_MU]   = m*c*c*pperp*pperp / (2.0*Babs);
        zval[GuidingCenterEquation::COORD_ZETA] = 0.0;

        p.SetCharge(q);
        p.SetMass(m);
        p.InitializeMomentum(
            Particle::COORDINATE_PPAR, Particle::COORDINATE_PPERP,
            ppar, pperp, 1.0, 1.0
        );
        p.InitializePosition(Particle::POSITION_GUIDINGCENTER, X0, 1.0, 0.0);

        eqn->InitializeParticle(&p, zval);
        eqn->Evaluate(0.0, zval, dzdt);

        Xdot[0] = dzdt[GuidingCenterEquation::COORD_X];
        Xdot[1] = dzdt[GuidingCenterEquation::COORD_Y];
        Xdot[2] = dzdt[GuidingCenterEquation::COORD_Z];

        return hypot(Xdot[0], Xdot[2]);
    };

    unsigned int i, n;
    slibreal_t a, b, c, d, h, yc, yd, L;
    const slibreal_t invphi = 0.5*(sqrt(5.0)-1.0), invphi2 = invphi*invphi;
    a = mf->GetMinRadius();
    b = mf->GetMaxRadius();
    h = b-a;
    L = mf->GetMagneticAxisR();

    n =(unsigned int)(ceil(log(tol/h) / log(invphi)));

    c = a + invphi2*h;
    d = a + invphi*h;

    yc = gxpol(c);
    yd = gxpol(d);

    for (i = 0; i < n; i++) {
        if (yc < yd) {
            b = d;
            d = c;
            yd = yc;
            h = invphi*h;
            c = a + invphi2*h;
            yc = gxpol(c);
        } else {
            a = c;
            c = d;
            yc = yd;
            h = invphi*h;
            d = a + invphi*h;
            yd = gxpol(d);
        }
    }

    if (yc < yd)
        return (0.5*(a+d) - L);
    else
        return (0.5*(c+b) - L);
}

/**
 * Check that the effective magnetic is found correctly.
 */
bool Test_ParticleGenerator::CheckEffectiveMagneticAxis() {
    unsigned int i, j;
    // Actual values are unimportant, except for the tolerance
    const string config =
        "a=0,1,2;\n"
        "ppar=10,12,2;\n"
        "pperp=1,2,2;\n"
        "driftshifttol="+to_string(10.0*REAL_EPSILON)+";\n";
    struct global_settings globset;
    globset.include_drifts = true;

    MagneticFieldAnalytical2D *mf = GetMagneticField();
    GuidingCenterEquation *eqn
        = new GuidingCenterEquation(
            mf, &globset
        );

    Configuration *conf = new Configuration();
    conf->FromString(config);
    ConfigBlock root = conf->GetRootBlock();
    ParticleGenerator *pg = new ParticleGenerator(mf, &root, &globset);
    Particle *p = pg->AllocateParticle();
    p->SetCharge(-ELEMENTARY_CHARGE);
    p->SetMass(ELECTRON_MASS);

    slibreal_t shift = 0.0, Delta = 0.0, rshift, ppar, pperp, dppar, dpperp;
    dppar  = (effax_ppar1  - effax_ppar0)  / (effax_nppar -1.0);
    dpperp = (effax_pperp1 - effax_pperp0) / (effax_npperp-1.0);
    for (i = 0; i < effax_nppar; i++) {
        ppar = effax_ppar0 + i*dppar;

        for (j = 0; j < effax_npperp; j++) {
            pperp = effax_pperp0 + i*dpperp;

            shift = pg->CalculateOrbitDriftShift(
                mf, p->GetMass(), p->GetCharge(), ppar, pperp
            );

            // Actual shift
            rshift = GetDriftShift(mf, eqn, ppar, pperp, pg->GetDriftShiftTolerance());

            Delta = fabs((rshift-shift) / rshift);
            if (isnan(Delta) || Delta > 10.0*pg->GetDriftShiftTolerance()) {
                this->PrintError("Effective magnetic axis with index [%u,%u] not calculated accurately. Delta = %e.", i, j, Delta);
                this->PrintError("ppar(%u) = %e, pperp(%u) = %e.", i, ppar, j, pperp);
                printf("rshift = %e, shift = %e\n", rshift, shift);

                delete mf;
                delete conf;
                delete pg;

                return false;
            }
        }
    }

    delete mf;
    delete conf;
    delete pg;

    return true;
}

/**
 * Check that the phase-space grid is properly
 * generated by the ParticleGenerator.
 */
bool Test_ParticleGenerator::CheckPhaseSpaceGrid() {
    slibreal_t a, b, c, d, e, f, t;
    const slibreal_t
        B0 = 5.0,
        rmajor = 0.68,
        rminor = 0.22,
        safety = 1.0;
    unsigned int i, j;
    MagneticFieldAnalytical2D *mf
        = new MagneticFieldAnalytical2D(B0, rmajor, rminor, MFATFS_CW, MFASF_CONSTANT, safety, 0.0);
    struct global_settings *globset = new struct global_settings;
    globset->include_drifts = true;

    const unsigned int NRCOORDS=3, NCOORDPAIRS=7;
    const char RCOORDS[NRCOORDS][4] =
        {"a","r","rho"};
    const slibreal_t RCOORDS_MIN[NRCOORDS] =
        {0.0,0.0,rmajor},
                     RCOORDS_MAX[NRCOORDS] =
        {1.0,rminor,rmajor+rminor};

    const char COORDPAIRS1[NCOORDPAIRS][6] =
        {"gamma","gamma","p","p","ppar","pperp","pperp"};
    const slibreal_t COORDPAIRS1_MIN[NCOORDPAIRS] =
        {1.0,    1.0,    0.0,0.0,0.0,   0.0,     0.0},
                     COORDPAIRS1_MAX[NCOORDPAIRS] =
        {100.0,  100.0,  100,100,100,   100.0,   100.0};
    const char COORDPAIRS2[NCOORDPAIRS][7] =
        {"thetap","xi","thetap","xi","pperp","thetap","xi"};
    const slibreal_t COORDPAIRS2_MIN[NCOORDPAIRS] =
        {0.0,     -1.0,0.0,     -1.0,0.0,    0.0,     -1.0},
                     COORDPAIRS2_MAX[NCOORDPAIRS] =
        {M_PI,    1.0,M_PI,     1.0, 100,    M_PI,    1.0};

    bool retval = true;
    for (i = 0; i < NRCOORDS; i++) {
        a = RCOORDS_MIN[i] + Rand() * (RCOORDS_MAX[i]-RCOORDS_MIN[i]);
        b = RCOORDS_MIN[i] + Rand() * (RCOORDS_MAX[i]-RCOORDS_MIN[i]);
        if (b < a) {t=b;b=a;a=t;}

        for (j = 0; j < NCOORDPAIRS; j++) {
            c = COORDPAIRS1_MIN[j] + Rand() * (COORDPAIRS1_MAX[j]-COORDPAIRS1_MIN[j]);
            d = COORDPAIRS1_MIN[j] + Rand() * (COORDPAIRS1_MAX[j]-COORDPAIRS1_MIN[j]);
            if (d < c) {t=d;d=c;c=t;}

            e = COORDPAIRS2_MIN[j] + Rand() * (COORDPAIRS2_MAX[j]-COORDPAIRS2_MIN[j]);
            f = COORDPAIRS2_MIN[j] + Rand() * (COORDPAIRS2_MAX[j]-COORDPAIRS2_MIN[j]);
            if (f < e) {t=f;f=e;e=t;}

            retval = Generate(
                mf, globset,
                RCOORDS[i], COORDPAIRS1[j], COORDPAIRS2[j],
                a, c, e, b, d, f
            );

            if (!retval) {
                delete globset;
                delete mf;
                return false;
            }

            // Test nr = np1 = 1
            retval = Generate(
                mf, globset,
                RCOORDS[i], COORDPAIRS1[j], COORDPAIRS2[j],
                a, c, e, a, c, e
            );

            if (!retval) {
                delete globset;
                delete mf;
                return false;
            }
        }
    }

    // ppar, pperp
    delete globset;
    delete mf;

    return true;
}

/**
 * Generates a phase-space grid according to the given
 * parameters and verifies that the 'ParticleGenerator'
 * class works.
 */
bool Test_ParticleGenerator::Generate(
    MagneticFieldAnalytical2D *mf, struct global_settings *globset,
    const string& r, const string& p1, const string& p2,
    const slibreal_t r0, const slibreal_t p10, const slibreal_t p20,
    const slibreal_t r1, const slibreal_t p11, const slibreal_t p21
) {
    bool retval = true;
    string config;
    char buffer1[50], buffer2[50];
    unsigned int i, j, k, nr, np1, np2;
    slibreal_t rarr[TEST_PARTICLE_GENERATOR_NRAD],
               p1arr[TEST_PARTICLE_GENERATOR_NP1],
               p2arr[TEST_PARTICLE_GENERATOR_NP2],
               Deltar, Deltap1, Deltap2,
               sr, sp1, sp2;

    nr = TEST_PARTICLE_GENERATOR_NRAD;
    np1 = TEST_PARTICLE_GENERATOR_NP1;
    np2 = TEST_PARTICLE_GENERATOR_NP2;

    if (r0 == r1) nr = 1;
    if (p10 == p11) np1 = 1;
    if (p20 == p21) np2 = 1;

    sprintf(buffer1, "%.16e", r0);
    sprintf(buffer2, "%.16e", r1);
    config = r+"="+string(buffer1)+","+string(buffer2)+","+to_string(nr)+";\n";

    sprintf(buffer1, "%.16e", p10);
    sprintf(buffer2, "%.16e", p11);
    config += p1+"="+string(buffer1)+","+string(buffer2)+","+to_string(np1)+";\n";

    sprintf(buffer1, "%.16e", p20);
    sprintf(buffer2, "%.16e", p21);
    config += p2+"="+string(buffer1)+","+string(buffer2)+","+to_string(np2)+";\n";

    /* Generate arrays */
    if (nr > 1) {
        for (i = 0; i < nr; i++)
            rarr[i] = r0 + (r1-r0) * ((slibreal_t)i)/((slibreal_t)(nr-1.0));
    } else rarr[0] = r0;

    if (np1 > 1) {
        for (i = 0; i < np1; i++)
            p1arr[i] = p10 + (p11-p10) * ((slibreal_t)i)/((slibreal_t)(np1-1.0));
    } else p1arr[0] = p10;

    if (np2 > 1) {
        for (i = 0; i < np2; i++)
            p2arr[i] = p20 + (p21-p20) * ((slibreal_t)i)/((slibreal_t)(np2-1.0));
    } else p2arr[0] = p20;

    Configuration *conf = new Configuration();
    conf->FromString(config);

    if (conf->HasError()) {
        this->PrintError("Configuration error: in text.\n%s", config.c_str());
        return false;
    }

    ConfigBlock root = conf->GetRootBlock();
    ParticleGenerator *pg = new ParticleGenerator(mf, &root, globset);
    Particle *p = pg->AllocateParticle();

    /* Generate particles */
    for (i = 0; i < np2; i++) {
        for (j = 0; j < np1; j++) {
            for (k = 0; k < nr; k++) {
                pg->Generate(p);

                sr  = GetValueFromParamName(p, r, mf);
                sp1 = GetValueFromParamName(p, p1, mf);
                sp2 = GetValueFromParamName(p, p2, mf);

                Deltar  = fabs((sr -rarr[k]) /rarr[k]);
                Deltap1 = fabs((sp1-p1arr[j])/p1arr[j]);
                Deltap2 = fabs((sp2-p2arr[i])/p2arr[i]);

                if (Deltap1 > 200.0*REAL_EPSILON) {
                    this->PrintError("ParticleGenerator: Generated momentum (%s) does not match expected value. Delta = %e.", p1.c_str(), Deltap1);
                    retval = false;
                } else if (Deltap2 > 200.0*REAL_EPSILON) {
                    this->PrintError("ParticleGenerator: Generated pitch angle (%s) does not match expected value. Delta = %e.", p2.c_str(), Deltap2);
                    retval = false;
                } else if (Deltar > 200.0*REAL_EPSILON) {
                    this->PrintError("ParticleGenerator: Generated radius (%s) does not match expected value. Delta = %e.", r.c_str(), Deltar);
                    retval = false;
                }

                if (!retval) {
                    delete conf;
                    delete p;
                    delete pg;
                    return false;
                }
            }
        }
    }

    if (!pg->IsFinished()) {
        this->PrintError("ParticleGenerator: All particles done, but particle generator not marked as 'finished'.");
        retval = false;
    }

    delete conf;
    delete p;
    delete pg;

    return retval;
}

/**
 * Get a parameter from the particle based on the
 * parameters name.
 */
slibreal_t Test_ParticleGenerator::GetValueFromParamName(Particle *p, const string& name, MagneticFieldAnalytical2D *mf) {
         if (name == "gamma")  return p->GetGamma();
    else if (name == "p")      return p->GetMomentum();
    else if (name == "ppar")   return p->GetPpar();
    else if (name == "pperp")  return p->GetPperp();
    else if (name == "thetap") return p->GetThetap();
    else if (name == "xi")     return p->GetXi();
    else if (name == "rho")    return p->GetRho();
    else if (name == "a") {
        slibreal_t r0 = mf->GetMagneticAxisR();
        slibreal_t rm = mf->GetMaxRadius() - r0;
        return ((p->GetRho()-r0) / rm);
    } else if (name == "r") {
        slibreal_t r0 = mf->GetMagneticAxisR();
        return (p->GetRho()-r0);
    } else
        throw SOFTLibException("GetValueFromParamName: Unrecognized parameter name: %s.", name.c_str());
}

/**
 * Check that invalid input is rejected.
 */
bool Test_ParticleGenerator::TestInvalidInput() {
    unsigned int i;
    const slibreal_t
        B0 = 5.0,
        rmajor = 0.68,
        rminor = 0.22,
        safety = 1.0;
    MagneticFieldAnalytical2D *mf
        = new MagneticFieldAnalytical2D(B0, rmajor, rminor, MFATFS_CW, MFASF_CONSTANT, safety, 0.0);
    struct global_settings *globset = new struct global_settings;
    globset->include_drifts = true;

    const unsigned int NA = 3;
    string valid_a = "a=0,1,2;";
    string invalid_a[NA] = {
        "a=0,1,0;",         // [0] Zero points
        "a=0,1,1;",         // [1] Different r0 and r1, but nr = 1
    };

    const unsigned int NP = 4;
    string valid_p = "p=1,10,1;";
    string invalid_p[NP] = {
        "p=1,10,0;",        // [0] Zero points
        "p=-0.1,10,2;",     // [1] Negative momentum (sign is set in pitch angle)
        "p=1,10,1;",        // [2] Different r0 and r1, but nr = 1
        "p=0,1,2;"          // [3] Zero momentum
    };

    string valid_thetap = "thetap=0.1,0.2,2;";
    string cnfstring;

    // Test invalid a
    Configuration *conf;
    for (i = 0; i < NA; i++) {
        cnfstring = invalid_a[i] + "\n" + valid_p + "\n" + valid_thetap;

        try {
            conf = new Configuration();
            conf->FromString(cnfstring);

            ConfigBlock root = conf->GetRootBlock();
            ParticleGenerator pg(mf, &root, globset);

            this->PrintError("ParticleGenerator invalid input #%u for 'a' did NOT throw an exception.", i);
            return false;
        } catch (SOFTLibException& ex) { }
    }

    // Test invalid p
    for (i = 0; i < NP; i++) {
        cnfstring = valid_a + "\n" + invalid_p[i] + "\n" + valid_thetap;

        try {
            conf = new Configuration();
            conf->FromString(cnfstring);

            ConfigBlock root = conf->GetRootBlock();
            ParticleGenerator pg(mf, &root, globset);

            this->PrintError("ParticleGenerator invalid input #%u for 'p' did NOT throw an exception.", i);
            return false;
        } catch (SOFTLibException& ex) { }
    }

    return true;
}

/**
 * Run the unit test for 'ParticleGenerator'.
 */
bool Test_ParticleGenerator::Run(bool) {
    bool success = true;
    if (CheckPhaseSpaceGrid())
        this->PrintOK("Phase-space grid generated properly.");
    else success = false;

    if (CheckEffectiveMagneticAxis())
        this->PrintOK("Effective magnetic axis calculated accurately.");
    else success = false;

    if (TestInvalidInput())
        this->PrintOK("Invalid input rejected.");
    else success = false;

    return success;
}

