#!/usr/bin/env python3
#
# ####################

import sys
if __name__ == '__main__':
    sys.path.append('../..')

import broadutil
import h5py
import numpy as np
import os.path
import time


TOLERANCE = 1e-7

# Quadrature resolution parameters
NR   = 400
NZ   = 400
NPHI = 400

# Maximum r/a to integrate to
MAXRA = 0.95


def run(soft):
    """
    Run this test

    soft: Path to SOFT executable
    """
    global TOLERANCE

    # Run SOFT
    print('Running SOFT...')
    broadutil.runSOFT("<pi>")

    # Load all output
    img = broadutil.Image("isotropic.mat")
    imgval = img.image[0,0]
    p = img.param1[0,0]
    sinThetap = np.sin(img.param2[0,0])

    # Compute received radiation using Python
    print('Computing integral...')

    t = time.time()
    recrad = ComputeReceivedRadiation(img)
    duration = time.time() - t

    # Add momentum jacobian (2pi from gyro integral)
    recrad *= p**2 * sinThetap / (2.0 * np.pi)

    print('Took {0}s'.format(duration))

    # Copmute relative error
    diso = np.abs((imgval - recrad)/imgval)

    success = True
    if diso > TOLERANCE:
        print('imgval = {0}'.format(imgval))
        print('recrad = {0}'.format(recrad))
        print('ERROR: SOFT does not integrate isotropic radiation properly. Delta = {0}'.format(diso))
        success = False

    return success


def ComputeReceivedRadiation(img):
    """
    Compute the total amount of radiation received by
    the detector used to generate 'img' (parameters should
    be within the object), if all particles (distributed
    uniformly in the plasma) emit radiation isotropically.
    """
    global NR, NZ, NPHI, MAXRA

    nhat  = img.detectorDirection[0,:]
    ehat1 = img.detectorEhat1[0,:]
    ehat2 = img.detectorEhat2[0,:]
    wall  = img.wall

    X0 = img.detectorPosition[0,:]
    l  = img.detectorAperture[0,0]

    A   = lambda l1, l2, cn, c1, c2: np.arctan((c1-l1)*(c2-l2)/(cn*np.sqrt(cn**2 + (c1-l1)**2 + (c2-l2)**2)))
    App = lambda cn, c1, c2: A(+l, +l, cn, c1, c2)
    Apn = lambda cn, c1, c2: A(+l, -l, cn, c1, c2)
    Anp = lambda cn, c1, c2: A(-l, +l, cn, c1, c2)
    Ann = lambda cn, c1, c2: A(-l, -l, cn, c1, c2)

    # Get domain boundaries
    rmin, rmax = np.amin(wall[0,:]), np.amax(wall[0,:])
    zmin, zmax = np.amin(wall[1,:]), np.amax(wall[1,:])

    minorRadius = 0.5*(rmax-rmin)
    majorRadius = 0.5*(rmax+rmin)
    radlim = MAXRA*MAXRA * minorRadius**2

    R   = np.linspace(rmin, rmax, NR)
    Z   = np.linspace(zmin, zmax, NZ)
    PHI = np.linspace(0, 2*np.pi, NPHI+1)
    PHI = PHI[:-1]

    ONE = np.ones(PHI.shape)
    X0 = np.array([x*ONE for x in X0])

    dR   = R[1]-R[0]
    dZ   = Z[1]-Z[0]
    dPHI = PHI[1]-PHI[0]

    COSPHI = np.cos(PHI)
    SINPHI = np.sin(PHI)

    I = 0.0
    for r in R:
        for z in Z:
            # Only include points that are within the domain
            mr = (r-majorRadius)**2 + z**2
            if mr > radlim:
                continue

            # Evaluate PHI integral as vectors
            x = np.array([r*COSPHI, r*SINPHI, z*ONE])
            rcp = x - X0

            cn = np.dot(nhat,  rcp)
            c1 = np.dot(ehat1, rcp)
            c2 = np.dot(ehat2, rcp)

            I += np.sum(App(cn, c1, c2) + Ann(cn, c1, c2) - Apn(cn, c1, c2) - Anp(cn, c1, c2))

    I = I * dR * dZ * dPHI

    return I


if __name__ == '__main__':
    if len(sys.argv) == 1:
        broadutil.init()
        sys.exit(0 if run('../../../../build/src/soft') else -1)
    else:
        print('ERROR: Unrecognized number of arguments.')
        sys.exit(-1)

