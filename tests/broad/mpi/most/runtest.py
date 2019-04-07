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

IMAGE_TOLERANCE = 1e-9
SPECTRUM_TOLERANCE = 1e-9
TOPVIEW_TOLERANCE = 1e-9

def run(soft):
    """
    Run this test

    soft: Path to SOFT executable
    """
    global IMAGE_TOLERANCE, SPECTRUM_TOLERANCE, TOPVIEW_TOLERANCE

    # Run with MPI
    print('Running SOFT with MPI...')
    broadutil.runSOFT_mpi("<pi>\nnum_threads=1;")

    # Load all output
    img = broadutil.Image("image.mat")
    imgval_mpi = np.sum(img.image)

    spc = broadutil.Spectrum("spectrum.mat")
    spcval_mpi = np.sum(spc.I)

    tov = broadutil.Image("topview.mat")
    tovval_mpi = np.sum(tov.image)

    # Run without MPI
    print('Running SOFT without MPI...')
    broadutil.runSOFT("<pi>")

    # Load all output
    img = broadutil.Image("image.mat")
    imgval = np.sum(img.image)

    spc = broadutil.Spectrum("spectrum.mat")
    spcval = np.sum(spc.I)

    tov = broadutil.Image("topview.mat")
    tovval = np.sum(tov.image)

    dimg = np.abs((imgval - imgval_mpi)/imgval)
    dspc = np.abs((spcval - spcval_mpi)/spcval)
    dtov = np.abs((tovval - tovval_mpi)/tovval)

    success = True
    if dimg > IMAGE_TOLERANCE:
        print('ERROR: Image was not computed properly when using MPI.')
        success = False
    if dspc > SPECTRUM_TOLERANCE:
        print('ERROR: Spectrum was not computed properly when using MPI.')
        success = False
    if dtov > TOPVIEW_TOLERANCE:
        print('ERROR: Topview was not computed properly when using MPI.')
        success = False

    return success


if __name__ == '__main__':
    if len(sys.argv) == 1:
        broadutil.init()
        sys.exit(0 if run('../../../build/src/soft') else -1)
    else:
        print('ERROR: Unrecognized number of arguments.')
        sys.exit(-1)
