/**
 * Initialization of the magnetic field.
 */

#include <vector>

#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include <softlib/MagneticField/MagneticFieldNumeric2D.h>

#include "Init/Init.h"
#include "Init/InitMagneticField.h"
#include "SOFT.h"
#include "SOFTException.h"

using namespace std;

/**
 * Initialize a magnetic field object from the
 * given ConfigBlock object.
 *
 * conf: Configuration of magnetic field.
 */
MagneticField2D *InitMagneticField(ConfigBlock *conf) {
    MagneticField2D *mf;
    Setting *type;

    if (!conf->HasSetting("type"))
        throw SOFTException("Type of magnetic field '%s' not specified.", conf->GetName().c_str());
    
    type = conf->GetSetting("type");
    if (type->GetNumberOfValues() != 1)
        throw SOFTException("Magnetic field '%s': type: Invalid value assigned to parameter. Expected single string.", conf->GetName().c_str());

    if (type->GetString() == "analytical")
        mf = InitMagneticFieldAnalytical(conf);
    else if (type->GetString() == "numeric")
        mf = InitMagneticFieldNumeric(conf);
    else
        throw SOFTException("Magnetic field '%s': type: Unrecognized magnetic field type: %s.", conf->GetName().c_str(), type->GetString().c_str());

    return mf;
}

/**
 * Initialize an analytical magnetic field based
 * on the given configuration.
 *
 * conf: Configuration of magnetic field.
 */
MagneticFieldAnalytical2D *InitMagneticFieldAnalytical(ConfigBlock *conf) {
    MagneticFieldAnalytical2D *mf;
    Setting *s;
    slibreal_t B0, Rm, rminor, qa1=1.0, qa2=1.0;
    enum MFAToroidalFieldSign tfs;
    enum MFASafetyFactorType qtype;
    string parent = "Magnetic field '"+conf->GetName()+"'";

    B0     = init_get_scalar(conf, "B0", parent);
    Rm     = init_get_scalar(conf, "Rm", parent);
    rminor = init_get_scalar(conf, "rminor", parent);

    // toroidal field sign
    if (!conf->HasSetting("sigmaB"))
        tfs = MFATFS_CW;
    else {
        if ((*conf)["sigmaB"] == "cw" || (*conf)["sigmaB"] == "+")
            tfs = MFATFS_CW;
        else if ((*conf)["sigmaB"] == "ccw" || (*conf)["sigmaB"] == "-")
            tfs = MFATFS_CCW;
        else
            throw SOFTException("%s: sigmaB: Invalid value assigned to parameter.", parent.c_str());
    }

    // qtype
    if (!conf->HasSetting("qtype"))
        throw SOFTException("%s: Required parameter 'qtype' not specified.", parent.c_str());
    s = conf->GetSetting("qtype");
    if (s->GetNumberOfValues() != 1)
        throw SOFTException("%s: qtype: Invalid value assigned to parameter. Expected single value.", parent.c_str());

    if (s->GetString() == "constant")
        qtype = MFASF_CONSTANT;
    else if (s->GetString() == "exponential")
        qtype = MFASF_EXPONENTIAL;
    else if (s->GetString() == "linear")
        qtype = MFASF_LINEAR;
    else if (s->GetString() == "quadratic")
        qtype = MFASF_QUADRATIC;
    else
        throw SOFTException("%s: qtype: Unrecognized q-profile specified: %s.", parent.c_str(), s->GetString().c_str());

    // q-parameters (optional)
    if (conf->HasSetting("qa1"))
        qa1 = init_get_scalar(conf, "qa1", parent);
    if (conf->HasSetting("qa2"))
        qa2 = init_get_scalar(conf, "qa2", parent);

    mf = new MagneticFieldAnalytical2D(B0, Rm, rminor, tfs, qtype, qa1, qa2, conf->GetName(), conf->GetName());
    return mf;
}

/**
 * Initialize an analytical magnetic field based
 * on the given configuration.
 *
 * conf: Configuration of magnetic field.
 */
MagneticFieldNumeric2D *InitMagneticFieldNumeric(ConfigBlock *conf) {
    MagneticFieldNumeric2D *mf;
    enum sfile_type ftype = SFILE_TYPE_UNKNOWN;
    string filename, filetype;
    string parent = "Magnetic field '"+conf->GetName()+"'";

    filename = init_get_string(conf, "filename", parent);

    // File type
    if (conf->HasSetting("filetype")) {
        filetype = init_get_string(conf, "filetype", parent);
        ftype = SFile::GetFileType(filetype);
    }

    mf = new MagneticFieldNumeric2D(filename, ftype);
    return mf;
}

