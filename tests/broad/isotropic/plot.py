#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt

h5f = h5py.File('isotropic.mat')

IMG = np.transpose(h5f['image'][:,:])

plt.imshow(IMG)
plt.show()
