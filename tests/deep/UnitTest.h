#ifndef _UNITTEST_H
#define _UNITTEST_H

#include <string>
#include <softlib/config.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>

using namespace std;

class UnitTest {
	protected:
		string name;

        const slibreal_t
            MF_B0 = 5.0,
            MF_Rm = 0.68,
            MF_rminor = 0.22;
	public:
		UnitTest(const string&);
		string& GetName();
		bool HasName(const string&);

        MagneticFieldAnalytical2D *GetMagneticField();

		void PrintError(const string&, ...);
		void PrintOK(const string&, ...);
		void PrintStatus(const string&, ...);
		void PrintWarning(const string&, ...);
		slibreal_t Rand();
		virtual bool Run(bool) = 0;
        void SeedRand();
};

#endif/*_UNITTEST_H*/
