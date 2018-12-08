#ifndef _RADIATION_SPACE3D_H
#define _RADIATION_SPACE3D_H

#include <string>
#include <softlib/config.h>

namespace __Radiation {
    class Space3D : public RadiationOutput {
        private:
            std::string output;
            Vector<3> point0, point1;
            long long int pixelsX, pixelsY, pixelsZ;

            slibreal_t *s3dimage=nullptr;
            size_t imagesize;

        public:
            Space3D(Detector *d, MagneticField2D *m) : RadiationOutput(d,m) {}

            void GetPoint(ConfigBlock*, const std::string&, Vector<3>&);
            void Initialize();
    };

    class Space3DException : public RadiationOutputException {
        public:
            template<typename ... Args>
            Space3DException(const std::string &msg, Args&& ... args)
                : RadiationOutputException(msg, std::forward<Args>(args) ...) {
                AddModule("Space3D");
            }
    };
}

#endif/*_RADIATION_SPACE3D_H*/
