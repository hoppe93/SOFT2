/**
 * Radiation :: Output :: SoVVolume
 *
 * Configuration of the 'SoVVolume' radiation output module.
 */

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include "Tools/Radiation/Output/SoVVolume.h"

using namespace __Radiation;
using namespace std;

const string SoVVolume::DEFAULT_QUANTITIES[] = {"none"};
/*template<typename T, unsigned int sz>
unsigned int __def_size(T(&)[sz]) { return sz; }
const unsigned int SoVVolume::NDEFAULT_QUANTITIES = __def_size(SoVVolume::DEFAULT_QUANTITIES);*/
const unsigned int SoVVolume::NDEFAULT_QUANTITIES = 0;

/**
 * Allocate the array holding the SoV volume.
 */
void SoVVolume::AllocateVolume() {
    this->volumearray = new slibreal_t[this->arraysize];
    for (size_t i = 0; i < this->narrayelements; i++)
        this->volumearray[i] = 0.0;
}

/**
 * Configure the 'SoVVolume' radiation output module.
 * 
 * conf: Configuration block for the module.
 * root: Root configuration block, providing access
 *       to the ConfigBlock's of possible sub-modules (unused).
 */
void SoVVolume::Configure(ConfigBlock *conf, ConfigBlock *__UNUSED__(root)) {
    this->SetName(conf->GetName());

	// common
	if (conf->HasSetting("common"))
		this->ConfigureCommonQuantities(DEFAULT_QUANTITIES, NDEFAULT_QUANTITIES, conf->GetSetting("common")->GetTextVector());
	else
		this->ConfigureCommonQuantities(DEFAULT_QUANTITIES, NDEFAULT_QUANTITIES);

    // output
    if (!conf->HasSetting("output"))
        throw SoVVolumeException("No output name has been specifed. Set the option 'output'.");
    else
        this->output = (*conf)["output"];

    this->narrayelements = this->np1*this->np2;
    this->arraysize      = this->narrayelements * sizeof(slibreal_t);
}

/**
 * Initialize the 3D image.
 */
void SoVVolume::Initialize() {
    AllocateVolume();
}

