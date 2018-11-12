/**
 * Print a welcome message for the 'Radiation' tool.
 */

#include <string>
#include "SOFT.h"
#include "Tools/Radiation/Radiation.h"

using namespace std;
using namespace __Radiation;

void Radiation::Welcome(const string &prefix) {
    // Detector
    SOFT::PrintInfo(prefix+"Detector:        "+this->detector->GetName().c_str());
    this->detector->
          PrintInfo(prefix+"                 ");

    // Model
    SOFT::PrintInfo(prefix+"Model:           "+this->model->GetDescription());
    SOFT::PrintInfo(prefix+"Polarization:    "+(this->MeasuresPolarization()?SOFT::PRINT_YES:SOFT::PRINT_NO));

    // Outputs
    if (noutput > 0)
        SOFT::PrintInfo(prefix+"Enabled outputs:");
    else
        SOFT::PrintWarning(SOFT::WARNING_TR_NO_OUTPUT, "No output modules have been coupled.");

    for (unsigned int i = 0; i < noutput; i++) {
        SOFT::PrintInfo(prefix+prefix+this->output[i]->GetName());
    }
}
