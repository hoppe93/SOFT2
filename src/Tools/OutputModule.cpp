/**
 * Implementation of routines common to all output modules.
 */

#include <iostream>
#include <string>
#include <vector>
#include "SOFT.h"
#include "Tools/OutputModule.h"

//#include <H5Cpp.h>

using namespace std;

const char
    OutputModule::DATETIME[]           = "datetime",
	OutputModule::RO_DOMAIN[]		   = "domain",
	OutputModule::DISTFUNC[]		   = "distribution",
	OutputModule::PARAM1[]             = "param1",
	OutputModule::PARAM1NAME[]         = "param1name",
	OutputModule::PARAM2[]             = "param2",
	OutputModule::PARAM2NAME[]         = "param2name",
	OutputModule::R[]                  = "r",
    OutputModule::RO_SEPARATRIX[]      = "separatrix",
	OutputModule::RO_WALL[]            = "wall",
    OutputModule::VERSION_SOFT[]       = "versionSOFT";
    //OutputModule::VERSION_SOFTLIB[]    = "versionSOFTLib";

/**
 * Define a common quantity.
 */
void OutputModule::DefineCommonQuantity(
    const string &name, const qfunc &func
) {
	#define PAIR std::pair<string, qfunc>
    all_quantities.insert(PAIR(name, func));
}

/**
 * Initialize the map containing all possible common
 * quantities that can be exported.
 */
void OutputModule::InitializeCommonQuantities() {

    // datetime
    DefineCommonQuantity(
        DATETIME, [this](SFile *sf) {
            struct tm *timeinfo;
            time_t rawtime;

            time(&rawtime);
            timeinfo = localtime(&rawtime);

            const int buffer_size = 200;
            char buffer[buffer_size];
            strftime(buffer, buffer_size, "%Y-%m-%d %H:%M:%S", timeinfo);

            sf->WriteString(this->DATETIME, buffer);
        }
    );

	// distribution
	DefineCommonQuantity(
		DISTFUNC, [this](SFile *sf) {
			slibreal_t *r = this->particlegenerator->GetRGrid();
			slibreal_t *param1 = this->particlegenerator->GetP1Grid();
			slibreal_t *param2 = this->particlegenerator->GetP2Grid();
			int p1type = this->particlegenerator->GetP1Type();
			int p2type = this->particlegenerator->GetP2Type();

			size_t nr = this->particlegenerator->GetNr();
			size_t n1 = this->particlegenerator->GetN1();
			size_t n2 = this->particlegenerator->GetN2();

			Particle *p = new Particle();
			slibreal_t **rhoeff = this->particlegenerator->GetRhoEff();

			DistributionFunction *dist = this->soft->distribution;
			slibreal_t *f = new slibreal_t[nr*n1*n2];
			for (size_t k = 0; k < nr; k++) {
				 for (size_t j = 0; j < n2; j++) {
				 	for (size_t i = 0; i < n1; i++) {
						p->InitializeMomentum(
							(enum Particle::coordinate)p1type,
							(enum Particle::coordinate)p2type,
							param1[i], param2[j],
							// zeta, dp1, dp2, dzeta
							0, 0, 0, 0
						);

						slibreal_t d = 0;
						if (this->particlegenerator->HasDrifts())
							d = rhoeff[i][j] - this->particlegenerator->GetRhoMin();

						f[(k*n2 + j)*n1 + i] = dist->Eval(r[k], p->GetMomentum(), p->GetXi(), d);
					}
				 }
			}

			sfilesize_t dims[3] = {nr, n2, n1};
			sf->WriteMultiArray("distribution", f, 3, dims);

			delete p;
			delete [] f;
		}
	);

	// param1
    DefineCommonQuantity(
		PARAM1, [this](SFile *sf) {
			slibreal_t *param1 = this->particlegenerator->GetP1Grid();
			unsigned int n1 = this->particlegenerator->GetN1();

			sf->WriteList(this->PARAM1, param1, n1);
		}
	);

	// param1name
    DefineCommonQuantity(
		PARAM1NAME, [this](SFile *sf) {
			int param1type = this->particlegenerator->GetP1Type();
			string param1name = Particle::GetCoordinateName(param1type);

			sf->WriteString("param1name", param1name);
		}
	);

	// param2
    DefineCommonQuantity(
		PARAM2, [this](SFile *sf) {
			slibreal_t *param2 = this->particlegenerator->GetP2Grid();
			unsigned int n2 = this->particlegenerator->GetN2();

			sf->WriteList(this->PARAM2, param2, n2);
		}
	);

	// param2name
    DefineCommonQuantity(
		PARAM2NAME, [this](SFile *sf) {
			int param2type = this->particlegenerator->GetP2Type();
			string param2name = Particle::GetCoordinateName(param2type);

			sf->WriteString("param2name", param2name);
		}
	);

	// r
    DefineCommonQuantity(
		R, [this](SFile *sf) {
			slibreal_t *r = this->particlegenerator->GetRGrid();
			unsigned int nr = this->particlegenerator->GetNr();

			sf->WriteList(this->R, r, nr);
		}
	);

	// domain, separatrix & wall
	DefineCommonQuantity(RO_DOMAIN, [this](SFile *sf) { this->write_domain(sf, this->RO_DOMAIN); });
    DefineCommonQuantity(RO_SEPARATRIX, [this](SFile *sf) { this->write_separatrix(sf, this->RO_SEPARATRIX); });
	DefineCommonQuantity(RO_WALL,   [this](SFile *sf) { this->write_domain(sf, this->RO_WALL); });

    // version
    DefineCommonQuantity(VERSION_SOFT, [this](SFile *sf) { sf->WriteString(this->VERSION_SOFT, SOFT_GIT_SHA1); });
    //DefineCommonQuantity(VERSION_SOFTLIB, [this](SFile *sf) { sf->WriteString(this->VERSION_SOFTLIB, SOFTLIB_GIT_SHA1); });
}

/**
 * Write the domain used to the given SFile.
 * Give the variable the given name 'name'.
 *
 * sf:   SFile object to write domain to.
 * name: Name of variable to store domain under.
 */
void OutputModule::write_domain(SFile *sf, const string &name) {
	unsigned int ndom = this->magfield->GetNDomain();
	slibreal_t **domain = new slibreal_t*[2], *rdom, *zdom;
	domain[0] = new slibreal_t[2*ndom];
	domain[1] = domain[0] + ndom;

	rdom = this->magfield->GetRDomain();
	zdom = this->magfield->GetZDomain();

	for (unsigned int i = 0; i < ndom; i++) {
		domain[0][i] = rdom[i];
		domain[1][i] = zdom[i];
	}

	sf->WriteArray(name, domain, 2, ndom);
}

/**
 * Write the separatrix to the given SFile
 * (if the magnetic field has a separatrix).
 *
 * sf:   SFile object to write separatrix to.
 * name: Name of variable to store domain under.
 */
void OutputModule::write_separatrix(SFile *sf, const string &name) {
	unsigned int nsep = this->magfield->GetNSeparatrix();
	slibreal_t **sep = new slibreal_t*[2], *rsep, *zsep;
	sep[0] = new slibreal_t[2*nsep];
	sep[1] = sep[0] + nsep;

	rsep = this->magfield->GetRSeparatrix();
	zsep = this->magfield->GetZSeparatrix();

	for (unsigned int i = 0; i < nsep; i++) {
		sep[0][i] = rsep[i];
		sep[1][i] = zsep[i];
	}

	sf->WriteArray(name, sep, 2, nsep);
}

/**
 * Interpret the given list of strings and build the
 * list of quantities to write to file.
 *
 * qList:     User-specified list of quantities to write.
 * defaults:  Default quantities for this OutputModule module.
 * ndefaults: Number of strings in the 'defaults' list.
 */
void OutputModule::ConfigureCommonQuantities(
	const string *defaults, const unsigned int ndefaults,
	const vector<string> &qList
) {
	common_quantities.clear();

	// If user-specified list is empty, treat as 'default'.
	if (qList.empty()) {
		for (unsigned int i = 0; i < ndefaults; i++)
			common_quantities.push_back(defaults[i]);
	}

	vector<string>::iterator
		cqb = common_quantities.begin(),
		cqe = common_quantities.end();

	for (vector<string>::const_iterator it = qList.begin(); it != qList.end(); it++) {
		// Append all available options
		if (*it == "all") {
			for (map<string, qfunc>::iterator it = all_quantities.begin(); it != all_quantities.end(); it++)
				common_quantities.push_back(it->first);
		// Append the defaults
		} else if (*it == "default") {
			for (unsigned int i = 0; i < ndefaults; i++)
				common_quantities.push_back(defaults[i]);
		// Clear list
		} else if (*it == "none") {
			common_quantities.clear();
		} else {
			// Remove quantities
			if (it->front() == '-') {
				vector<string>::iterator el = cqb;
				string item = it->substr(1);

				// Locate element
				while (el != cqe && *el != item) el++;

				// Remove
				if (el != cqe)
					common_quantities.erase(el);
			// Append quantities
			} else {
				string s;
				if (it->front() == '+')
					s = it->substr(1);
				else
					s = *it;

				// Make sure element doesn't already exist
				vector<string>::iterator el = cqb;
				while (el != cqe && *el != s) el++;

				if (el == cqe)
					common_quantities.push_back(s);
			}
		}
	}
}

/**
 * Write all common quantities named in 'qList' to the
 * output file represented by 'output'.
 *
 * output: Output 'SFile' to which quantities should be written.
 */
void OutputModule::WriteCommonQuantities(SFile *output) {
	vector<string> notFound;

		for (vector<string>::iterator it = common_quantities.begin(); it != common_quantities.end(); it++) {
			if (all_quantities.find(*it) != all_quantities.end()) {
				all_quantities[*it](output);
			} else
				notFound.push_back(*it);
		}

	// Emit warning if some quantity was not found
	if (notFound.size() > 0) {
		string acum = "";
		for (vector<string>::iterator it = notFound.begin(); it != notFound.end(); it++) {
			if (it != notFound.begin())
				acum += ", ";

			acum += *it;
		}

		SOFT::PrintWarning(SOFT::WARNING_TRO_COMMON_NOT_RECOGNIZED, "%s: The following common quantities were not recognized: %s.", this->name.c_str(), acum.c_str());
	}
}

