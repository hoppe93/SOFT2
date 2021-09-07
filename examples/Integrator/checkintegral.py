#!/usr/bin/env python3
#
# This script runs SOFT on the 'pi' file located in the same directory as this
# script, parses the output and fetches the value of the resulting integral.
# This value is then compared to the analytically known value of the volume for
# a torus (to which the unit distribution function should integrate).
#
# NOTE: This script must be run from the directory in which the file is located.

import numpy as np
import subprocess
import sys


def main():
    p = subprocess.Popen(['../../build/src/soft', 'pi'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outp = p.communicate()

    stdout_data = outp[0].decode('utf-8')
    stderr_data = outp[1].decode('utf-8')

    if p.returncode != 0:
        print(stderr_data)
        raise RuntimeError('SOFT exited with a non-zero exit code.')

    V_soft = 0
    lines = stdout_data.split('\n')
    for line in lines:
        if line.strip()[:9] == 'Integral:':
            V_soft = float(line.split(' ')[-1])
            break

    # Magnetic field parameters
    R = 0.68        # Major radius
    a = 0.22 * 0.5  # Minor radius
    V = 2*np.pi**2 * R * a**2

    print('V_soft   = {},  V = {}'.format(V_soft, V))
    print('Ratio    = {}'.format(V_soft/V))

    return 0


if __name__ == '__main__':
    sys.exit(main())

