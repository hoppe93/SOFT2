/**
 * Implementation of routines common to all output modules.
 */

#include <string>
#include <vector>
#include "Tools/Radiation/Output/RadiationOutput.h"

using namespace __Radiation;
using namespace std;

const char
	RadiationOutput::DETECTOR_APERTURE[]  = "detectorAperture",
	RadiationOutput::DETECTOR_DIRECTION[] = "detectorDirection",
	RadiationOutput::DETECTOR_POSITION[]  = "detectorPosition",
	RadiationOutput::DETECTOR_VISANG[]    = "detectorVisang",
	RadiationOutput::RO_DOMAIN[]		  = "domain",
	RadiationOutput::RO_WALL[]            = "wall",
	RadiationOutput::TP_BOUNDARY[]		  = "tpBoundary";

/**
 * Initialize the map containing all possible common
 * quantities that can be exported.
 */
void RadiationOutput::InitializeCommonQuantities() {
	#define PAIR std::pair<string, qfunc>
	// detectorAperture
	all_quantities.insert(
		PAIR(DETECTOR_APERTURE, [this](SFile *sf) {
			slibreal_t detap = this->detector->GetAperture();
			sf->WriteList(this->DETECTOR_APERTURE, &detap, 1);
		})
	);

	// detectorDirection
	all_quantities.insert(
		PAIR(DETECTOR_DIRECTION, [this](SFile *sf) {
			slibreal_t detdir[3];
			this->detector->GetDirection().ToArray(detdir);
			sf->WriteList(this->DETECTOR_DIRECTION, detdir, 3);
		})
	);

	// detectorPosition
	all_quantities.insert(
		PAIR(DETECTOR_POSITION, [this](SFile *sf) {
			slibreal_t detpos[3];
			this->detector->GetPosition().ToArray(detpos);
			sf->WriteList(this->DETECTOR_POSITION, detpos, 3);
		})
	);

	// detectorVisang
	all_quantities.insert(
		PAIR(DETECTOR_VISANG, [this](SFile *sf) {
			slibreal_t visang = 2.0 * this->detector->GetVisionAngleFOV();
			sf->WriteList(this->DETECTOR_VISANG, &visang, 1);
		})
	);

	// domain & wall
	all_quantities.insert(PAIR(RO_DOMAIN, [this](SFile *sf) { this->write_domain(sf, this->RO_DOMAIN); }));
	all_quantities.insert(PAIR(RO_WALL,   [this](SFile *sf) { this->write_domain(sf, this->RO_WALL); }));
}

/**
 * Write the domain used to the given SFile.
 * Give the variable the given name 'name'.
 *
 * sf:   SFile object to write domain to.
 * name: Name of variable to store domain under.
 */
void RadiationOutput::write_domain(SFile *sf, const string &name) {
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
 * Interpret the given list of strings and build the
 * list of quantities to write to file.
 *
 * qList:     User-specified list of quantities to write.
 * defaults:  Default quantities for this RadiationOutput module.
 * ndefaults: Number of strings in the 'defaults' list.
 */
void RadiationOutput::ConfigureCommonQuantities(
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

				if (el != cqe)
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
void RadiationOutput::WriteCommonQuantities(SFile *output) {
	for (vector<string>::iterator it = common_quantities.begin(); it != common_quantities.end(); it++) {
		if (all_quantities.find(*it) != all_quantities.end())
			all_quantities[*it](output);
	}
}

