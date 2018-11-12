# Utility functions for the broad testing script

import os.path
import sys
import subprocess
import h5py

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

class Image:
    def __init__(self, filename):
        with h5py.File(filename, 'r') as f:
            self.image = f['image'][:,:]
            self.detectorPosition = f['detectorPosition'][:,:]
            self.detectorDirection = f['detectorDirection'][:,:]
            self.detectorVisang = f['detectorVisang'][:,:]
            self.wall = f['wall'][:,:]

