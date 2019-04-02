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
TMIN = 0.1
TMAX = 0.3

def loadGreensFunction(filename):
    GF = None
    with h5py.File(filename, "r") as f:
        func = f['func'][:,:]
        frmt = "".join(map(chr, f['type'][:,:][:,0].tolist()))

        if frmt != "12ij":
            raise Exception("Invalid format of Greens's function: "+frmt);

        #GF = np.reshape(func, [2,2,300,300])
        GF = func.shape

        param1name = "".join(map(chr, f['param1name'][:,:][:,0].tolist()))
        param2name = "".join(map(chr, f['param2name'][:,:][:,0].tolist()))

        p, thetap = None, None
        if param1name == 'p' and param2name == 'thetap':
            p = f['param1'][:,:]
            thetap = f['param2'][:,:]
        elif param1name == 'thetap' and param1name == 'p':
            p = f['param2'][:,:]
            thetap = f['param1'][:,:]
        else:
            raise Exception("Green's function coordinates were not p and thetap.")

        P, THETAP = np.meshgrid(p, thetap)
        diffel = P * P * np.sin(THETAP)
    
    return GF, diffel
    
def plotGreen(i,j):
    """
    Plot image (i,j) of Green's function
    """
    GF = loadGreensFunction("green.mat")

    plt.contourf(GF[i,j,:,:])
    plt.show()

def runGreen():
    global PMIN, PMAX, TMIN, TMAX
    broadutil.runSOFT("<pi>\n@Radiation rad { output = green; }\n@ParticleGenerator PGen { p = "+str(PMIN)+","+str(PMAX)+",2; thetap = "+str(TMIN)+","+str(TMAX)+",2; }")

    # Load Green's function
    green, diffel = loadGreensFunction("green.mat")
    return green, diffel

def runImage(p, thetap):
    broadutil.runSOFT("<pi>\n@Radiation rad { output = image; }\n@ParticleGenerator PGen { p = "+str(p)+","+str(p)+",1; thetap = "+str(thetap)+","+str(thetap)+",1; }")

    # Load image
    img = broadutil.Image("image.mat")

    return img
    
def run(soft):
    """
    Run this test

    soft: Path to SOFT executable
    """
    global PMIN, PMAX, TMIN, TMAX

    # Run Green's function
    print("Generating Green's function...")
    green, diffel = runGreen()

    # Run verification images
    i = 0
    P = [PMIN, PMAX]
    T = [TMIN, TMAX]
    for p in [0,1]:
        for t in [0,1]:
            i += 1
            print("Generating image #{0}...".format(i))
            img = runImage(P[p], T[t])
            I   = green[p,t,:,:] * diffel[p,t]
            if np.amax(np.abs(img.image-I)) > np.sqrt(np.finfo(float).eps):
                print("ERROR: Green's function element [{0},{1}] did not match the corresponding image.".format(p,t))
                return False

    return True


if __name__ == '__main__':
    if len(sys.argv) == 1:
        broadutil.init()
        run('../../../../build/src/soft')
    elif len(sys.argv) == 2:
        if sys.argv[1] == 'plot':
            plotGreen(0,0)

