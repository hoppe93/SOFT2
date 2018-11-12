#!/usr/bin/env python3
#
# Script for generating lookup tables of SOFT orbits.
#
# HOW TO RUN:
#   ./genorbits.py               -- Generate 60 GC orbits
#                                   (30 with drifts, 30 without)
#                                   & 30 particle orbits.
#   ./genorbits.py n             -- Generate 2n GC orbits
#                                   (n with drifts, n without)
#                                   & n particle orbits.
#   ./genorbits.py nGC np        -- Generate 2nGC GC orbits
#                                   (nGC with drifts, nGC without)
#                                   & np particle orbits.
#   ./genorbits.py nGC nGCnd np  -- Generate nGC GC orbits with
#                                   drifts, nGCnd GC orbits without
#                                   drifts, and np particle orbits.
# ###################################################################

import io, os, sys
import numpy as np
import scipy.io
import subprocess
import random

SOFTPATH = '.'
OUTPUT = 'orbits.cpp'

try:
    SOFTPATH = os.environ['SOFTPATH']

    # Make sure path doesn't end with slash
    if SOFTPATH[-1] == '/':
        SOFTPATH = SOFTPATH[:-1]
except KeyError:
    print('WARNING: Unable to determine SOFT path.')
    SOFTPATH = None

def runSOFT(pifile):
    """
    Run SOFT, passing the given pifile on stdin.
    The contents of stderr are interpreted as CSV output
    and are returned as an array.
    """
    global SOFTPATH

    if SOFTPATH is None:
        raise RuntimeError('The path to SOFT has not been specified.')
    
    p = subprocess.Popen([SOFTPATH+'/soft'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

    stderr_data = p.communicate(input=bytearray(pifile, 'ascii'))[1].decode('utf-8')

    if p.returncode is not 0:
        raise RuntimeError('SOFT exited with a non-zero exit code.')

    c = io.StringIO(stderr_data)
    return np.loadtxt(c, delimiter=',', skiprows=1)

def generatePi(r, p, thetap, particleOrbit=False, drifts=True):
    """
    Generate SOFTv1 pi file

    r:             Particle initial position.
    p:             Particle momentum.
    thetap:        Particle pitch angle.
    particleOrbit: If True, solves the guiding-center equations
                   of motion. Otherwise, the particle equations
                   of motion.
    """
    pi = ""
    c = 299792458.0
    time = 2*np.pi*r*np.sqrt(p**2 + 1) / (c*p*np.cos(thetap))

    if particleOrbit:
        pi += "useequation=particle-relativistic;\n"
    else:
        pi += "useequation=guiding-center-relativistic;\n"

    pi += "usetool=orbit;\n"
    if drifts: pi += "nodrifts=no;\n"
    else: pi += "nodrifts=yes;\n"

    pi += "tolerance=1e-12;\n"

    pi += "magnetic_field=circular;\n"
    pi += "magnetic circular { B0=5; major_radius=0.68; minor_radius=0.22; safety_factor=1; }\n"

    pi += "particles {\n"
    
    if particleOrbit: pi += "    gc_position=no;\n"
    else: pi += "    gc_position=yes;\n"

    pi += "    t=0,"+str(time)+";\n"
    pi += "    r="+str(r)+","+str(r)+",1;\n"
    pi += "    p="+str(p)+","+str(p)+",1;\n"
    pi += "    pitch="+str(thetap)+","+str(thetap)+",1;\n"
    pi += "}\n"

    pi += "tool orbit {\n"
    pi += "    output=@stderr;\n"
    pi += "}\n"

    return pi

def generateOrbit(particleOrbit, drifts):
    """
    Generate an orbit.
    """
    ELECTRON_MASS_EV = 0.5109989461 * 1e6

    # Limits
    RMIN, RMAX = 0.69, 0.90
    PMIN, PMAX = 1e6, 5e7
    TMIN, TMAX = 0, 1.4

    rr = random.random()
    rp = random.random()
    rt = random.random()

    r = RMIN + (RMAX-RMIN)*rr
    p = PMIN + (PMAX-PMIN)*rp
    thetap = TMIN + (TMAX-TMIN)*rt

    out = runSOFT(generatePi(r, p, thetap, particleOrbit, drifts))
    return {
        'r': r,
        'p': p / ELECTRON_MASS_EV,
        'thetap': thetap,
        't': out[:,0],
        'x': out[:,1],
        'y': out[:,2],
        'z': out[:,3]
    }


def generateOrbits(nGC, nGCnd, nP):
    """
    Generate guiding-center and particle orbits with
    SOFTv1 that can be used to test SOFTv2. A C-file
    for inclusion into a test program is generated.

    nGC:   Number of GC orbits (w/ drifts)
    nGCnd: Number of GC orbits (w/o drifts)
    nP:    Number of particles orbits
    """
    GCorbits = []
    GCNDorbits = []
    Porbits = []

    GCorbits_max = 0
    GCNDorbits_max = 0
    Porbits_max = 0

    print("GC orbit (w/ drifts):   ", end="")
    initProgressBar(nGC)
    for i in range(0, nGC):
        o = generateOrbit(False, True)
        GCorbits.append(o)
        GCorbits_max = max(GCorbits_max, len(o['t']))
        progress(i)

    print("GC orbit (w/o drifts):  ", end="")
    initProgressBar(nGCnd)
    for i in range(0, nGCnd):
        o = generateOrbit(False, False)
        GCNDorbits.append(o)
        GCNDorbits_max = max(GCNDorbits_max, len(o['t']))
        progress(i)

    print("Particle orbits:        ", end="")
    initProgressBar(nP)
    for i in range(0, nP):
        o = generateOrbit(True, True)
        Porbits.append(o)
        Porbits_max = max(Porbits_max, len(o['t']))
        progress(i)

    # Output
    generateOutput(GCorbits, GCNDorbits, Porbits, GCorbits_max, GCNDorbits_max, Porbits_max)

def appendOrbitData(orbit):
    frm = "{:.9e}"

    l = len(orbit['t'])
    s = 1
    if l > 100:
        s = int(np.floor(l / 100))

    n = len(range(0, l, s))

    cfile = "\n    {"
    cfile += frm.format(orbit['r'])+","
    cfile += frm.format(orbit['p'])+","
    cfile += frm.format(orbit['thetap'])+","

    arrt = "(slibreal_t["+str(n)+"]){"
    arrx = "(slibreal_t["+str(n)+"]){"
    arry = "(slibreal_t["+str(n)+"]){"
    arrz = "(slibreal_t["+str(n)+"]){"

    for i in range(0, l, s):
        arrt += frm.format(orbit['t'][i])+","
        arrx += frm.format(orbit['x'][i])+","
        arry += frm.format(orbit['y'][i])+","
        arrz += frm.format(orbit['z'][i])+","

    arrt = arrt[:-1]+"},"
    arrx = arrx[:-1]+"},"
    arry = arry[:-1]+"},"
    arrz = arrz[:-1]+"}"

    cfile += str(n)+","
    cfile += arrt+arrx+arry+arrz+"},"

    return cfile
    
def generateOutput(GCorbits, GCNDorbits, Porbits, GCmax, GCNDmax, Pmax):
    """
    Generate a C++ file with orbits.
    """
    global OUTPUT
    cfile = "/* This file was autogenerated by 'genorbits.py' */\n\n"

    cfile += "#include <softlib/config.h>\n"
    cfile += "#include \"Test_ParticlePusher.h\"\n\n"

    # GC orbits
    cfile += "const unsigned int Test_ParticlePusher::ngcorbits = "+str(len(GCorbits))+";\n"
    cfile += "const unsigned int Test_ParticlePusher::guiding_center_orbits_nmax = "+str(GCmax)+";\n"
    cfile += "struct orbit Test_ParticlePusher::guiding_center_orbits[Test_ParticlePusher::ngcorbits] = {"

    for orbit in GCorbits:
        cfile += appendOrbitData(orbit)

    cfile = cfile[:-1]
    cfile += "\n};\n\n"

    # GC orbits (w/o drifts)
    cfile += "const unsigned int Test_ParticlePusher::ngcorbits_nd = "+str(len(GCNDorbits))+";\n"
    cfile += "const unsigned int Test_ParticlePusher::guiding_center_orbits_nd_nmax = "+str(GCNDmax)+";\n"
    cfile += "struct orbit Test_ParticlePusher::guiding_center_orbits_nd[Test_ParticlePusher::ngcorbits_nd] = {"

    for orbit in GCNDorbits:
        cfile += appendOrbitData(orbit)

    cfile = cfile[:-1]
    cfile += "\n};\n\n"

    # Particle orbits
    cfile += "const unsigned int Test_ParticlePusher::npartorbits = "+str(len(Porbits))+";\n"
    cfile += "const unsigned int Test_ParticlePusher::particle_orbits_nmax = "+str(Pmax)+";\n"
    cfile += "struct orbit Test_ParticlePusher::particle_orbits[Test_ParticlePusher::npartorbits] = {"

    for orbit in Porbits:
        cfile += appendOrbitData(orbit)

    cfile = cfile[:-1]
    cfile += "\n};\n\n"

    # Write file
    with open(OUTPUT, 'w') as cppf:
        cppf.write(cfile)

    print("Wrote output to '"+OUTPUT+"'. Done.")

def initProgressBar(mx):
    global PROGRESSMAX, PROGRESSWIDTH, PROGRESSINDEX

    sys.stdout.write("[%s]" % ("-" * PROGRESSWIDTH))
    sys.stdout.flush()
    sys.stdout.write("\b" * (PROGRESSWIDTH+1))

    PROGRESSMAX = mx 
    PROGRESSINDEX = 0

def progress(v):
    global PROGRESSMAX, PROGRESSWIDTH, PROGRESSINDEX

    if PROGRESSINDEX == PROGRESSWIDTH: return

    w = int(np.round(PROGRESSWIDTH * ((v+1)/PROGRESSMAX)))

    if PROGRESSINDEX >= w: return

    while PROGRESSINDEX < w:
        sys.stdout.write("#")
        PROGRESSINDEX += 1

    sys.stdout.flush()

    if PROGRESSINDEX == PROGRESSWIDTH:
        sys.stdout.write("]\n")
        sys.stdout.flush()

PROGRESSMAX = 100
PROGRESSWIDTH = 30
PROGRESSINDEX = 0

if __name__ == '__main__':
    nGC = 30
    nGCnd = 30
    nP  = 30

    if len(sys.argv) == 2:
        nGC = nGCnd = nP = int(sys.argv[1])
    elif len(sys.argv) == 3:
        nGC = nGCnd = int(sys.argv[1])
        nP = int(sys.argv[2])
    elif len(sys.argv) == 4:
        nGC = int(sys.argv[1])
        nGCnd = int(sys.argv[2])
        nP = int(sys.argv[3])

    random.seed()

    generateOrbits(nGC, nGCnd, nP)
    
