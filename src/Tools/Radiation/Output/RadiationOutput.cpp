/**
 * Implementation of common methods for '@RadiationOutput' modules.
 */

#include <iostream>
#include <string>
#include <vector>
#include "SOFT.h"
#include "Tools/OutputModule.h"
#include "Tools/Radiation/Output/RadiationOutput.h"

using namespace __Radiation;
using namespace std;


const char 
	RadiationOutput::DETECTOR_APERTURE[]  = "detectorAperture",
	RadiationOutput::DETECTOR_DIRECTION[] = "detectorDirection",
	RadiationOutput::DETECTOR_EHAT1[]     = "detectorEhat1",
	RadiationOutput::DETECTOR_EHAT2[]     = "detectorEhat2",
	RadiationOutput::DETECTOR_POSITION[]  = "detectorPosition",
    RadiationOutput::DETECTOR_ROLL[]      = "detectorRoll",
	RadiationOutput::DETECTOR_VISANG[]    = "detectorVisang";

void RadiationOutput::InitializeCommonQuantities() {
    // Initialize base common quantities first
    this->OutputModule::InitializeCommonQuantities();

	// detectorAperture
    DefineCommonQuantity(
		DETECTOR_APERTURE, [this](SFile *sf) {
			slibreal_t detap = this->detector->GetAperture();
			sf->WriteList(this->DETECTOR_APERTURE, &detap, 1);
		}
	);

	// detectorDirection
    DefineCommonQuantity(
		DETECTOR_DIRECTION, [this](SFile *sf) {
			slibreal_t detdir[3];
			this->detector->GetDirection().ToArray(detdir);
			sf->WriteList(this->DETECTOR_DIRECTION, detdir, 3);
		}
	);

	// detectorEhat1
    DefineCommonQuantity(
		DETECTOR_EHAT1, [this](SFile *sf) {
			slibreal_t e1[3];
			this->detector->GetEHat1().ToArray(e1);
			sf->WriteList(this->DETECTOR_EHAT1, e1, 3);
		}
	);

	// detectorEhat2
    DefineCommonQuantity(
		DETECTOR_EHAT2, [this](SFile *sf) {
			slibreal_t e2[3];
			this->detector->GetEHat2().ToArray(e2);
			sf->WriteList(this->DETECTOR_EHAT2, e2, 3);
		}
	);

	// detectorPosition
    DefineCommonQuantity(
		DETECTOR_POSITION, [this](SFile *sf) {
			slibreal_t detpos[3];
			this->detector->GetPosition().ToArray(detpos);
			sf->WriteList(this->DETECTOR_POSITION, detpos, 3);
		}
	);

	// detectorRoll
    DefineCommonQuantity(
        DETECTOR_ROLL, [this](SFile *sf) {
            slibreal_t roll = this->detector->GetRoll();
            sf->WriteList(this->DETECTOR_ROLL, &roll, 1);
        }
    );

	// detectorVisang
    DefineCommonQuantity(
		DETECTOR_VISANG, [this](SFile *sf) {
			slibreal_t visang = 2.0 * this->detector->GetVisionAngleFOV();
			sf->WriteList(this->DETECTOR_VISANG, &visang, 1);
		}
	);
}

