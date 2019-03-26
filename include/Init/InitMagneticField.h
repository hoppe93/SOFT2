#ifndef _INIT_MAGNETIC_FIELD_H
#define _INIT_MAGNETIC_FIELD_H

#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>
#include <softlib/MagneticField/MagneticFieldAnalytical2D.h>
#include <softlib/MagneticField/MagneticFieldLUKE.h>
#include <softlib/MagneticField/MagneticFieldNumeric2D.h>

MagneticField2D *InitMagneticField(ConfigBlock*);
MagneticFieldAnalytical2D *InitMagneticFieldAnalytical(ConfigBlock*);
MagneticFieldLUKE *InitMagneticFieldLUKE(ConfigBlock*);
MagneticFieldNumeric2D *InitMagneticFieldNumeric(ConfigBlock*);

#endif/*_INIT_MAGNETIC_FIELD_H*/
