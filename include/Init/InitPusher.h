#ifndef _INIT_PUSHER_H
#define _INIT_PUSHER_H

#include <string>
#include <softlib/config.h>
#include <softlib/Configuration.h>
#include <softlib/MagneticField/MagneticField2D.h>

#include "Orbit/ParticlePusher.h"
#include "SOFT.h"

ParticlePusher *InitPusher(MagneticField2D*, struct global_settings*, Configuration*);
string InitPusherDefaults();

#endif/*_INIT_PUSHER_H*/
