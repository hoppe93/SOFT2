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


TOLERANCE = 1e-7

# Quadrature resolution parameters
NR   = 400
NZ   = 400
NPHI = 100


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

    # Compute received radiation using Python
    print('Computing integral...')
    recrad = ComputeReceivedRadiation(img)

    # Copmute relative error
    diso = np.abs((imgval - recrad)/imgval)

    print('imgval = {0}'.format(imgval))
    print('recrad = {0}'.format(recrad))

    success = True
    if diso > TOLERANCE:
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
    global NR, NZ

    nhat  = img.detectorDirection
    ehat1 = img.detectorEhat1
    ehat2 = img.detectorEhat2
    wall  = img.wall

    X0 = img.detectorPosition
    l  = img.detectorAperture

    App = lambda cn, c1, c2: np.arctan((c1-l)*(c2-l)/(cn*np.sqrt(cn**2 + (c1-l)**2 + (c2-l)**2)))
    Anp = lambda cn, c1, c2: np.arctan((c1+l)*(c2-l)/(cn*np.sqrt(cn**2 + (c1+l)**2 + (c2-l)**2)))
    Apn = lambda cn, c1, c2: np.arctan((c1-l)*(c2+l)/(cn*np.sqrt(cn**2 + (c1-l)**2 + (c2+l)**2)))
    Ann = lambda cn, c1, c2: np.arctan((c1+l)*(c2+l)/(cn*np.sqrt(cn**2 + (c1+l)**2 + (c2+l)**2)))

    # Get domain boundaries
    rmin, rmax = np.amin(wall[0,:]), np.amax(wall[0,:])
    zmin, zmax = np.amin(wall[1,:]), np.amax(wall[1,:])

    minorRadius = 0.5*(rmax-rmin)
    minorRadius2 = minorRadius**2

    R   = np.linspace(rmin, rmax, NR)
    Z   = np.linspace(zmin, zmax, NZ)
    PHI = np.linspace(0, 2*np.pi, NPHI+1)
    PHI = PHI[:-1]

    dR   = R[1]-R[0]
    dZ   = Z[1]-Z[0]
    dPHI = PHI[1]-PHI[0]

    COSPHI = np.cos(PHI)
    SINPHI = np.sin(PHI)

    I = 0.0
    for r in R:
        for z in Z:
            # Only include points that are within the domain
            mr = r**2 + z**2
            if mr > minorRadius2:
                continue

            for i in range(0, len(PHI)):
                x = np.array([r*COSPHI[i], r*SINPHI[i], z])
                rcp = X0 - x

                cn = np.dot(nhat,  rcp)
                c1 = np.dot(ehat1, rcp)
                c2 = np.dot(ehat2, rcp)

                I += App(cn, c1, c2) + Ann(cn, c1, c2) - Apn(cn, c1, c2) - Anp(cn, c1, c2)
                
    I = dR * dZ * dPHI

    return I


if __name__ == '__main__':
    if len(sys.argv) == 1:
        broadutil.init()
        sys.exit(0 if run('../../../../build/src/soft') else -1)
    else:
        print('ERROR: Unrecognized number of arguments.')
        sys.exit(-1)
