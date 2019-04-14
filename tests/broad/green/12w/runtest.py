#!/usr/bin/env python3
#
# Verifies that Green's functions of
# the format '12ij' contain the same
# images as are generated with just
# the 'image' sycout.
# #####################################

import sys
if __name__ == '__main__':
    sys.path.append('../..')

import broadutil
import h5py
import numpy as np
import numpy.matlib
import os.path
import matplotlib.pyplot as plt

PMIN = 20.0
PMAX = 40.0
NP   = 2
TMIN = 0.1
TMAX = 0.3
NT   = 2
TOLERANCE = 1e-7

def plotGreen(i,j):
    """
    Plot image (i,j) of Green's function
    """
    GF = broadutil.Green("green.mat")

    plt.plot(GF.func[i,j,:])
    plt.show()

def runGreen():
    global PMIN, PMAX, NP, TMIN, TMAX, NT
    broadutil.runSOFT("<pi>\n@Radiation rad { output = green; }\n@ParticleGenerator PGen { p = "+str(PMIN)+","+str(PMAX)+","+str(NP)+"; thetap = "+str(TMIN)+","+str(TMAX)+","+str(NT)+"; }")

    # Load Green's function
    green = broadutil.Green("green.mat")
    P, THETAP = np.meshgrid(green.param1, green.param2)
    diffel = P**2 * np.sin(THETAP)
    return green, diffel

def runSpectrum(p, thetap):
    broadutil.runSOFT("<pi>\n@Radiation rad { output = spectrum; }\n@ParticleGenerator PGen { p = "+str(p)+","+str(p)+",1; thetap = "+str(thetap)+","+str(thetap)+",1; }")

    # Load spectrum
    return broadutil.Spectrum("spectrum.mat")
    
def run(soft):
    """
    Run this test

    soft: Path to SOFT executable
    """
    global PMIN, PMAX, NP, TMIN, TMAX, NT, TOLERANCE

    # Run Green's function
    print("Generating Green's function...")
    green, diffel = runGreen()

    # Run verification spectra
    i = 0
    P = np.linspace(PMIN, PMAX, NP)
    T = np.linspace(TMIN, TMAX, NT)
    for p in range(0, NP):
        for t in range(0, NT):
            i += 1
            print("Generating spectrum #{0}...".format(i))
            spec = runSpectrum(P[p], T[t])
            I    = green.func[p,t,:] * diffel[p,t]

            dspec = np.amax(np.abs((spec.I-I)/I))
            if dspec > TOLERANCE:
                print("ERROR: Green's function element [{0},{1}] did not match the corresponding spectrum. Delta = {2}".format(p,t,dspec))
                return False

    return True


if __name__ == '__main__':
    if len(sys.argv) == 1:
        broadutil.init()
        sys.exit(run('../../../../build/src/soft'))
    elif len(sys.argv) == 3:
        i, j = int(sys.argv[1]), int(sys.argv[2])
        plotGreen(i, j)
    else:
        print('ERROR: Expected at most one command-line argument.')

