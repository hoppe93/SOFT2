#!/usr/bin/env python3
#
# Verify that the distribution functions are
# applied properly in SOFT.
#
# To generate a test distribution function, run
#
#   ./runtest.py gen
#
# NOTE: This test does not verify that all
# distribution functions are implemented correctly,
# but only that they are actually applied.
# ####################

import sys
if __name__ == '__main__':
    sys.path.append('..')

import broadutil
import h5py
import numpy as np
import numpy.matlib
import os.path

DISTRIBUTION = 'distribution.h5'
RMIN = 1.7
RMAX = 2.2
PMAX = 100.0

_a = np.array([5.0])
_p = np.array([50.0])
_rho = np.array([1.9])
_xi = np.array([.98])

def distributionFunction(a, r, p, xi):
    global RMIN, RMAX, PMAX

    if np.any(xi < -1.0) or np.any(xi > 1.0):
        raise RuntimeError('Invalid xi value encountered.')

    fr = (RMAX-r)/(RMAX-RMIN)
    fp = (PMAX-p)/PMAX
    fx = np.abs(xi)

    fr[np.where(r <  RMIN)] = 0
    fr[np.where(r >= RMAX)] = 0
    fp[np.where(p <  0)]    = 0
    fp[np.where(p >= PMAX)] = 0

    return a*fr*fp*fx

def generateDistribution():
    global RMIN, RMAX, PMAX, _a, _rho, _p, _xi, DISTRIBUTION

    a = _a
    NP  = 11
    NXI = 10
    NR  = 10
    p  = np.array([np.linspace(1, PMAX, NP)])
    xi = np.array([np.linspace(-1.0, 1.0, NXI)])
    r  = np.array([np.linspace(RMIN, RMAX, NR)]).T

    P, XI = np.meshgrid(p, xi)
    P  = np.matlib.repmat(np.reshape(P, (1,NP*NXI)), NR, 1)
    XI = np.matlib.repmat(np.reshape(XI, (1,NP*NXI)), NR, 1)
    
    R = np.matlib.repmat(r, 1, NP*NXI)

    F = distributionFunction(a, R, P, XI)

    fp0  = np.array([distributionFunction(a, np.array([r[0,0]]), p[0,:], np.array([xi[0,0]]))])
    fxi0 = np.array([distributionFunction(a, np.array([r[0,0]]), np.array([p[0,0]]), xi[0,:])])
    fr0  = np.array([distributionFunction(a, r[:,0], np.array([p[0,0]]), np.array([xi[0,0]]))])

    with h5py.File(DISTRIBUTION, 'w') as f:
        f.create_dataset('r', data=r.T)
        f.create_dataset('p', data=p)
        f.create_dataset('xi', data=xi)

        f.create_dataset('f', data=F)
        f.create_dataset('fp0', data=fp0)
        f.create_dataset('fxi0', data=fxi0)
        f.create_dataset('fr0', data=fr0)

        f.create_dataset('punits', data='normalized')

def run(soft):
    """
    Run this test

    soft: Path to SOFT executable
    """
    global _a, _rho, _p, _xi, DISTRIBUTION

    if not os.path.isfile(DISTRIBUTION):
        generateDistribution()

    broadutil.runSOFT("<pi>")
    img = broadutil.Image("image.mat")
    nval = np.sum(img.image)

    broadutil.runSOFT("<pi>\ndistribution_function = dist;")
    img = broadutil.Image("image.mat")
    ival = np.sum(img.image)

    f = distributionFunction(_a, _rho, _p, _xi)[0]
    d = np.abs(ival/nval - f)
    if d > 2*np.finfo(float).eps:
        print('ERROR: Distribution function was not applied properly. Delta = '+str(d))
        return False

    return True


if __name__ == '__main__':
    if len(sys.argv) == 1:
        broadutil.init()
        run('../../../build/src/soft')
    elif len(sys.argv) == 2:
        if sys.argv[1] == 'gen':
            generateDistribution()
        elif sys.argv[1] == 'eval':
            print(distributionFunction(_a, _rho, _p, _xi)[0])
            #print(distributionFunction(_a, np.array([1.75556]), np.array([1.0]), np.array([-1.0]))[0])

