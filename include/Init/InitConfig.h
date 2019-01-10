#ifndef _INIT_CONFIG_H
#define _INIT_CONFIG_H

#include <softlib/Configuration.h>

/* Confblock type IDs */
extern int
    CONFBLOCK_EQUATION,
    CONFBLOCK_MAGNETICFIELD,
    CONFBLOCK_PARTICLEGENERATOR,
    CONFBLOCK_PARTICLEPUSHER,
    CONFBLOCK_DISTRIBUTION,
    CONFBLOCK_RADIALPROFILE;

void InitConfig(Configuration*);

#endif/*_INIT_CONFIG_H*/
