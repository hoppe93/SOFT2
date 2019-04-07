#!/usr/bin/env python3
#
# Verify that SOFT generates the same output when
# run with MPI as when run without MPI.
# ####################

import sys
if __name__ == '__main__':
    sys.path.append('../..')

import broadutil
import h5py
import numpy as np
import os.path

TOLERANCE = 5e-7

def run(soft):
    """
    Run this test

    soft: Path to SOFT executable
    """
    global TOLERANCE

    # Run with MPI
    print('Running SOFT with MPI...')
    broadutil.runSOFT_mpi("<pi>\nnum_threads=1;")

    # Load output
    vol = broadutil.SoVVolume("sovvolume.mat")
    volval_mpi = np.sum(vol.volumearray)

    # Run without MPI
    print('Running SOFT without MPI...')
    broadutil.runSOFT("<pi>")

    # Load all output
    vol = broadutil.SoVVolume("sovvolume.mat")
    volval = np.sum(vol.volumearray)

    dvol = np.abs((volval - volval_mpi)/volval)

    success = True
    if dvol > TOLERANCE:
        print('ERROR: SoV volume was not computed properly when using MPI. Delta = {0:.5e}'.format(dvol))
        success = False

    return success


if __name__ == '__main__':
    if len(sys.argv) == 1:
        broadutil.init()
        sys.exit(0 if run('../../../build/src/soft') else -1)
    else:
        print('ERROR: Unrecognized number of arguments.')
        sys.exit(-1)
