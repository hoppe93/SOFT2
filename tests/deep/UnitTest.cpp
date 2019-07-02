/**
 * Implementation of UnitTest class.
 */

#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <string>

#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include "UnitTest.h"

using namespace std;

/**
 * Constructor
 */
UnitTest::UnitTest(const string& name) {
    this->name = name;
    SeedRand();
}

/**
 * Get name of this unit test.
 */
string& UnitTest::GetName() { return this->name; }

/**
 * Check if this unit test has
 * the given name.
 *
 * cmp: Name to check against.
 */
bool UnitTest::HasName(const string& cmp) { return (this->name==cmp); }

/********************************************
 * FUNCTIONS FOR CONSTRUCTING DUMMY OBJECTS *
 *                                          *
 * These objects can be used if you "just   *
 * need a random object" for testing and    *
 * don't care about the specific settings   *
 * of the object.                           *
 ********************************************/
/**
 * Constructs a new MagneticFieldAnalytical2D object.
 */
MagneticFieldAnalytical2D *UnitTest::GetMagneticField() {
    return new MagneticFieldAnalytical2D(
        MF_B0, MF_Rm, MF_zaxis, MF_rminor, MFAFS_CW, MFAFS_CCW, MFASF_CONSTANT, 1.0, 0.0
    );
}

/**
 * Print an [ERROR] message.
 * Uses printf syntax.
 */
void UnitTest::PrintError(const string& s, ...) {
    va_list args;
    va_start(args, s);

    fprintf(stderr, "\x1B[1;31m[ERROR]\x1B[0m   ");
    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}
/**
 * Print an [OK] message.
 * Uses printf syntax.
 */
void UnitTest::PrintOK(const string& s, ...) {
    va_list args;
    va_start(args, s);

    fprintf(stderr, "\x1B[1;32m[OK]\x1B[0m      --> ");
    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}
/**
 * Print a general status message.
 * Uses printf syntax.
 */
void UnitTest::PrintStatus(const string& s, ...) {
    va_list args;
    va_start(args, s);

    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}
/**
 * Print a [WARNING] message.
 * Uses printf syntax.
 */
void UnitTest::PrintWarning(const string& s, ...) {
    va_list args;
    va_start(args, s);

    fprintf(stderr, "\x1B[1;93m[WARNING]\x1B[0m ");
    vfprintf(stderr, s.c_str(), args);
    fprintf(stderr, "\n");

    va_end(args);
}

/**
 * Generate a random number between
 * 0.0 and 1.0 (inclusive in both limits).
 */
slibreal_t UnitTest::Rand() {
	return (((slibreal_t)rand()) / ((slibreal_t)RAND_MAX));
}

/**
 * Seed the random number generator.
 */
void UnitTest::SeedRand() {
    srand(time(NULL));
}

