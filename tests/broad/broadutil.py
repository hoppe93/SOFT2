# Utility functions for the broad testing script

import os.path
import sys
import subprocess
import h5py
import tempfile

SOFTPATH = '../../build/src/soft'

def error(msg):
    print('ERROR: '+msg)
    sys.exit(1)

def init():
    """
    Initialize this module.
    """
    global SOFTPATH

    if not os.path.isfile(SOFTPATH):
        SOFTPATH = '../'+SOFTPATH
        if not os.path.isfile(SOFTPATH):
            SOFTPATH = '../'+SOFTPATH
            if not os.path.isfile(SOFTPATH):
                raise RuntimeError('Unable to find SOFT executable.')
    
def runSOFT(pi):
    """
    Run SOFT, passing on 'pi' on stdin to the executable.
    """
    global SOFTPATH

    p = subprocess.Popen([SOFTPATH], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

    stderr_data = p.communicate(input=bytearray(pi, 'ascii'))[1].decode('utf-8')

    if p.returncode is not 0:
        raise RuntimeError("SOFT exited with a non-zero exit code ("+str(p.returncode)+").\n"+stderr_data)

    return stderr_data

def runSOFT_mpi(pi):
    """
    Run SOFT throught MPIRUN utility, passing
    on 'pi' on stdin to the executable.
    """
    global SOFTPATH

    pifile, path = tempfile.mkstemp()
    try:
        with os.fdopen(pifile, 'w') as f:
            f.write(pi)

        args = ["mpirun", "-n", "4", SOFTPATH, path]
        env = os.environ

        p = subprocess.Popen(args, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
        stderr_data = p.communicate()[0].decode('utf-8')
    finally:
        os.remove(path)

    if p.returncode is not 0:
        raise RuntimeError("SOFT exited with a non-zero exit code ("+str(p.returncode)+").\n{}".format(stderr_data))

    return

class SOFTOutputBase:
    def tostring(self, arr):
        """
        Convert the given MATLAB string to a Python string.
        """
        return "".join(map(chr, arr[:,:][:,0].tolist()))
    
class Image(SOFTOutputBase):
    def __init__(self, filename):
        with h5py.File(filename, 'r') as f:
            self.image = f['image'][:,:]
            self.detectorPosition = f['detectorPosition'][:,:]
            self.detectorDirection = f['detectorDirection'][:,:]
            self.detectorVisang = f['detectorVisang'][:,:]
            self.wall = f['wall'][:,:]

class Green(SOFTOutputBase):
    def __init__(self, filename):
        with h5py.File(filename, 'r') as f:
            self.func = f['func'][:]
            self.param1 = f['param1'][:]
            self.param2 = f['param2'][:]
            self.r      = f['r'][:]
            self.wavelengths = f['wavelengths'][:]

            self.param1name = self.tostring(f['param1name'])
            self.param2name = self.tostring(f['param2name'])
            self.format     = self.tostring(f['type'])

class SoVVolume(SOFTOutputBase):
    def __init__(self, filename):
        with h5py.File(filename, 'r') as f:
            self.volumearray = f['volumearray'][:]
            self.param1      = f['param1'][:]
            self.param2      = f['param2'][:]

            self.param1name  = self.tostring(f['param1name'])
            self.param2name  = self.tostring(f['param2name'])

class Spectrum(SOFTOutputBase):
    def __init__(self, filename):
        with h5py.File(filename, 'r') as f:
            self.I = f['I'][:,:]
            self.wavelengths = f['wavelengths'][:,:]

            if 'Q' in f: self.Q = f['Q'][:,:]
            else: self.Q = None

            if 'U' in f: self.U = f['U'][:,:]
            else: self.U = None

            if 'V' in f: self.V = f['V'][:,:]
            else: self.V = None
            
