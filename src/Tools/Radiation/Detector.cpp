/**
 * Implementation of a detector object.
 *
 * Specifies properties of the detector.
 */

#include <limits>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/Vector.h>
#include "Tools/Radiation/Detector.h"
#include "Tools/Radiation/DetectorException.h"
#include "Tools/Radiation/Radiation.h"

#include "Tools/Radiation/Optics/Korger.h"
#include "Tools/Radiation/Optics/Optics.h"

using namespace std;
using namespace __Radiation;

/**
 * Constructor.
 *
 * aperture:     Detector aperture (side length of square detector).
 * visang:       Detector vision angle (of field-of-view).
 * direction:    Detector viewing direction (not necessarily normalized).
 * position:     Detector position vector. Relative to tokamak point of symmetry.
 * nwavelengths: Number of wavelengths to measure in spectral range
 *               (0 = do not measure spectrum).
 * lambda0:      Lower wavelength limit in spectral range.
 * lambda1:      Upper wavelength limit in spectral range.
 */
__Radiation::Detector::Detector(
    slibreal_t aperture, slibreal_t visang, Vector<3>& direction,
    Vector<3>& position, unsigned int nwavelengths,
    slibreal_t lambda0, slibreal_t lambda1, Optics *optics
) {
    this->aperture = aperture;
    this->direction = direction;
    this->position = position;
    this->name = "<unknown>";
    this->SetVisionAngle(visang, VISANG_TYPE_FOV);
    this->SetSpectralRange(lambda0, lambda1, nwavelengths);
    this->optics = optics;

    this->direction.Normalize();

    SetupBasisVectors();
}

/**
 * Constructor.
 *
 * conf: Configuration block specifying the settings
 *       for this detector object.
 */
__Radiation::Detector::Detector(ConfigBlock *conf, ConfigBlock *root) {
    Setting *s;
    this->name = conf->GetName();

    // aperture
    if (!conf->HasSetting("aperture"))
        throw DetectorException("Detector '%s': No detector aperture has been specified.", this->name.c_str());
    s = conf->GetSetting("aperture");
    if (!s->IsScalar())
        throw DetectorException("Detector '%s': Invalid specification of detector aperture. Expected positive real number.", this->name.c_str());
    else {
        this->aperture = s->GetScalar();

        if (this->aperture <= 0)
            throw DetectorException("Detector '%s': Invalid value assigned to detector aperture. Expected positive real number.", this->name.c_str());
    }

    // direction
    if (!conf->HasSetting("direction"))
        throw DetectorException("Detector '%s': No detector viewing direction has been specified.", this->name.c_str());
    s = conf->GetSetting("direction");
    if (!s->IsNumericVector(3))
        throw DetectorException("Detector '%s': Invalid specification of vector viewing direction. Expected three (3) numeric values.", this->name.c_str());
    else {
        vector<slibreal_t> v = s->GetNumericVector();
        slibreal_t *a = v.data();
        this->direction = Vector<3>(a);

        if (this->direction.Norm() == 0)
            throw DetectorException("Detector '%s': Null-vector specified as detector viewing direction.", this->name.c_str());

        this->direction.Normalize();
    }

    // optics
    if (conf->HasSetting("optics")) {
        string oname = (*conf)["optics"];
        if (!root->HasSubBlock(Radiation::CONFBLOCK_T_DETECTOR_OPTICS, oname))
            throw DetectorException("Detector '%s': Unrecognized optics '%s' specified.", this->name.c_str(), s->GetString().c_str());

        SetOptics(root->GetConfigBlock(Radiation::CONFBLOCK_T_DETECTOR_OPTICS, oname));
    } else
        SetOptics(root->GetConfigBlock(Radiation::CONFBLOCK_T_DETECTOR_OPTICS, "__default__"));

    // position
    if (!conf->HasSetting("position"))
        throw DetectorException("Detector '%s': No detector position has been specified.", this->name.c_str());
    s = conf->GetSetting("position");
    if (!s->IsNumericVector(3))
        throw DetectorException("Detector '%s': Invalid specification of vector position. Expected three (3) numeric values.", this->name.c_str());
    else {
        vector<slibreal_t> v = s->GetNumericVector();
        slibreal_t *a = v.data();
        this->position = Vector<3>(a);
    }

    // vision angle
    if (!conf->HasSetting("vision_angle"))
        throw DetectorException("Detector '%s': No detector vision angle has been specified.", this->name.c_str());
    s = conf->GetSetting("vision_angle");
    if (s->GetNumberOfValues() == 1) {
        if (!s->IsScalar())
            throw DetectorException("Detector '%s': Invalid specification of vision angle: parameter is not a number.", this->name.c_str());
        this->SetVisionAngle(s->GetScalar(), VISANG_TYPE_FOV);
    } else if (s->GetNumberOfValues() == 2) {
        vector<string> args = s->GetTextVector();
        slibreal_t visang = stod(args.front(), NULL);
        string unit = args.back();

        if (unit == "fov")
            this->SetVisionAngle(visang, VISANG_TYPE_FOV);
        else if (unit == "image")
            this->SetVisionAngle(visang, VISANG_TYPE_IMAGE);
        else
            throw DetectorException("Detector '%s': Unrecognized unit for vision angle: %s.", this->name.c_str(), unit.c_str());
    } else
        throw DetectorException("Detector '%s': Invalid specification of vision angle: too many parameters.", this->name.c_str());

    // spectrum
    if (conf->HasSetting("spectrum")) {
        s = conf->GetSetting("spectrum");
        if (s->IsNumericVector(3)) {
            SetSpectralRange(
                s->GetScalar(0),
                s->GetScalar(1),
                s->GetUnsignedInteger32(2)
            );
        } else if (s->IsBool() && s->GetBool() == false) {
            SetSpectralRange(0.0, 0.0, 0);      // Disable spectral range by setting nwavelengths = 0
        } else
            throw DetectorException("Detector '%s': Invalid format of detector spectral range. Expected 'real,real,int' or 'no'.", this->name.c_str());
    }

    SetupBasisVectors();
}

/**
 * Configure the optics model to use.
 */
void __Radiation::Detector::SetOptics(ConfigBlock *conf) {
    if (conf->GetSecondaryType() == "korger")
        this->optics = new Korger(this, conf);
    else
        throw DetectorException("Detector '%s': Unrecognized type '%s' of optics model '%s'.", this->name.c_str(), conf->GetSecondaryType().c_str(), conf->GetName().c_str());
}

/**
 * Calculates r_cp, the vector from the detector to the
 * particle, located at 'x'.
 *
 * x: Particle position.
 * 
 * RETURNS r_cp = x-R, where R is the detector position.
 */
Vector<3> __Radiation::Detector::CalculateRCP(Vector<3>& x) {
    return (x - position);
}

/**
 * Prints information about this detector to stdout.
 * 
 * prefix: Prefix string to print.
 */
void __Radiation::Detector::PrintInfo(const std::string &prefix) const {
    SOFT::PrintInfo(
        prefix+"ehat1 = (%.3f, %.3f, %.3f)\n"+prefix+"ehat2 = (%.3f, %.3f, %.3f)",
        this->ehat1[0], this->ehat1[1], this->ehat1[2],
        this->ehat2[0], this->ehat2[1], this->ehat2[2]
    );
}

/**
 * Sets the spectral range and initializes the
 * detector wavelength vector.
 * 
 * w0: Lower limit of spectral range.
 * w1: Upper limit of spectral range.
 * nw: Number of points in the discretized spectral range.
 */
void __Radiation::Detector::SetSpectralRange(slibreal_t w0, slibreal_t w1, unsigned int nw) {
    this->wavelength0 = w0;
    this->wavelength1 = w1;
    this->nwavelengths = nw;

    if (nw > 1) {
        this->wavelengths = new slibreal_t[nw];
        for (unsigned int i = 0; i < nw; i++)
            this->wavelengths[i] = w0 + (w1-w0) * ((slibreal_t)i) / (nw-1.0);
    } else if (nw == 1) {
        if (w0 != w1)
            throw DetectorException("Detector '%s': Spectral range consists of one point but endpoints differ.", this->name.c_str());

        this->wavelengths = new slibreal_t;
        this->wavelengths[0] = w0;
    }
}

/**
 * Sets the vision angle of the camera. The
 * 'type' parameter which type of vision angle
 * that is specified. Either, the vision angle
 * for the field-of-view can be specified, or
 * one can specify the vision angle for the pyramid
 * inscribed in the field-of-view.
 * NOTE: The vision angle is interpreted as the
 * half-angle (not the whole angle as in SOFT 1).
 *
 * visang: Vision angle.
 * type:   Type of the given vision angle. Must be
 *         either 'Detector::VISANG_TYPE_FOV' or
 *         'Detector::VISANG_TYPE_IMAGE'.
 */
void __Radiation::Detector::SetVisionAngle(slibreal_t visang, int type) {
    switch (type) {
        case VISANG_TYPE_FOV:
            vision_angle_fov       = visang;
            cos_vision_angle_fov   = cos(visang);
            tan_vision_angle_fov   = tan(visang);
            tan_vision_angle_image = tan_vision_angle_fov / sqrt(2.0);
            vision_angle_image     = atan(tan_vision_angle_image);
            break;
        case VISANG_TYPE_IMAGE:
            vision_angle_image     = visang;
            tan_vision_angle_image = tan(visang);
            tan_vision_angle_fov   = sqrt(2.0) * tan_vision_angle_image;
            vision_angle_fov       = atan(tan_vision_angle_fov);
            cos_vision_angle_fov   = cos(vision_angle_fov);
            break;
        default:
            throw DetectorException("Detector '%s': Invalid type of vision angle: %d.\n", this->name.c_str(), type);
    }

    if (vision_angle_fov <= 0.0)
        throw DetectorException("Detector '%s': Vision angle must be positive.", this->name.c_str());
    else if (vision_angle_fov >= 0.5*M_PI)
        throw DetectorException("Detector '%s': Vision angle may not be greater than pi.", this->name.c_str());
}

/**
 * Calculate the basis vectors ehat1 & ehat2
 * which together with nhat ('direction') form
 * an orthonormal coordinate system that spans
 * the detector plane.
 */
void __Radiation::Detector::SetupBasisVectors() {
    /*if (direction[1] == 0) {
        ehat2[0] = ehat2[2] = 0;
        ehat2[1] = 1.0;
    } else {
        ehat2[0] = direction[1];
        ehat2[1] =-direction[0];

        ehat2.Normalize();
    }

    Vector<3>::Cross(ehat2, direction, ehat1);*/
    if (direction[1] == 0) {
        ehat1[0] = ehat1[2] = 0;
        ehat1[1] = 1.0;
    } else {
        ehat1[0] = direction[1];
        ehat1[1] =-direction[0];
        ehat1[2] = 0.0;

        ehat1.Normalize();
    }

    // Right-handed system
    Vector<3>::Cross(direction, ehat1, ehat2);
    /*ehat2[0] = direction[1]*ehat1[2] - direction[2]*ehat1[1];
    ehat2[1] = direction[2]*ehat1[0] - direction[0]*ehat1[2];
    ehat2[2] = direction[0]*ehat1[1] - direction[1]*ehat1[0];*/
    /*
    // Left-handed system (as in SOFTv1)
    ehat2[0] = direction[2]*ehat1[1] - direction[1]*ehat1[2];
    ehat2[1] = direction[0]*ehat1[2] - direction[2]*ehat1[0];
    ehat2[2] = direction[1]*ehat1[0] - direction[0]*ehat1[1];
    */
}

