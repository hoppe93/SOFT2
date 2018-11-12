#ifndef _DETECTOR_H
#define _DETECTOR_H

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/Vector.h>

namespace __Radiation {
    class Detector {
        private:
            slibreal_t aperture;
            Vector<3> direction;
            Vector<3> position;
            slibreal_t vision_angle_fov,
                vision_angle_image,
                cos_vision_angle_fov,
                tan_vision_angle_fov,
                tan_vision_angle_image;

            std::string name;

            slibreal_t wavelength0, wavelength1, *wavelengths;
            unsigned int nwavelengths=0;

            // Basis vectors
            Vector<3> ehat1, ehat2;

        public:
            Detector(slibreal_t, slibreal_t, Vector<3>&, Vector<3>&, unsigned int, slibreal_t l0=0, slibreal_t l1=0);
            Detector(ConfigBlock*);

            Vector<3> CalculateRCP(Vector<3>&);

            slibreal_t GetAperture() { return aperture; }
            Vector<3>& GetDirection() { return direction; }
            Vector<3>& GetEHat1() { return ehat1; }
            Vector<3>& GetEHat2() { return ehat2; }
            Vector<3>& GetPosition() { return position; }
            slibreal_t GetCosVisionAngleFOV() { return cos_vision_angle_fov; }
            slibreal_t GetVisionAngleFOV() { return vision_angle_fov; }
            slibreal_t GetVisionAngleImage() { return vision_angle_image; }
            slibreal_t GetTanVisionAngleFOV() { return tan_vision_angle_fov; }
            slibreal_t GetTanVisionAngleImage() { return tan_vision_angle_image; }
            slibreal_t GetWavelengthLower() { return wavelength0; }
            slibreal_t GetWavelengthUpper() { return wavelength1; }
            slibreal_t *GetWavelengths() { return wavelengths; }
            unsigned int GetNWavelengths() { return nwavelengths; }

            const std::string GetName() const { return this->name; }
            void PrintInfo(const std::string &prfx="") const;

            bool HasSpectralRange() { return (nwavelengths > 0); }
            void SetVisionAngle(slibreal_t, int);
            void SetupBasisVectors();
            void SetSpectralRange(slibreal_t, slibreal_t, unsigned int);

            // Specifying type of vision angle
            enum {
                VISANG_TYPE_FOV,            // Vision angle is for field-of-view.
                VISANG_TYPE_IMAGE           // Vision angle is for image.
            };
    };
}

#endif/*_DETECTOR_H*/
