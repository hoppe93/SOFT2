#!/usr/bin/env python3
#
# Plot the distribution function evaluated by SOFT.
#

import h5py
import matplotlib.pyplot as plt
import numpy as np
from DREAM import GeriMap


with h5py.File('12.h5', 'r') as f:
    dist = f['distribution'][:]
    p = f['param1'][:]
    xi = f['param2'][:]
    r = f['r'][:]

N = p.size+1
cmap = GeriMap.get(N=N)
for i in range(p.size):
    plt.semilogy(xi, dist[-1,:,i], color=cmap(i/N), label=f'$p = {p[i]}mc$')

plt.xlabel(r'$\xi$')
plt.ylabel(r'$f(\xi)$')
plt.legend()

plt.show()

