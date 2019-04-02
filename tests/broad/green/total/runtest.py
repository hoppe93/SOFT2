#!/usr/bin/env python3
#
# Verify that the unit momentum distribution with
# a linear radial profile is applied properly in SOFT.
# ####################

import sys
if __name__ == '__main__':
    sys.path.append('../..')

import broadutil
import h5py
import numpy as np
import numpy.matlib
import os.path

RMIN = 0
RMAX = 0.95
MINOR_RADIUS = 2.2

def run(soft):
    """
    Run this test

    soft: Path to SOFT executable
    """
    broadutil.runSOFT("<pi>")
    img = broadutil.Image("image.mat")
    nvalI = np.sum(img.image)

    gf = broadutil.Green("green.mat")
    dr = gf.r[0,1] - gf.r[0,0]
    nvalG = np.sum(gf.func) * dr

    print(dr)

    print(nvalI)
    print(nvalG)

    d = np.abs((nvalI - nvalG) / nvalI)
    if d > 2*np.finfo(float).eps:
        print('ERROR: Radial profile was not applied properly to unit distribution function. Delta = '+str(d))
        return False

    return True


if __name__ == '__main__':
    broadutil.init()
    succ = run('../../../../build/src/soft')

    if succ:
        sys.exit(0)
    else:
        sys.exit(1)

