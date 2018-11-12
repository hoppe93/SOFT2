#!/usr/bin/env python3

import numpy as np
import h5py
import matplotlib.pyplot as plt

h5f = h5py.File('cone.mat')

IMG = np.transpose(h5f['image'][:,:])

if np.count_nonzero(IMG) == 0:
    print('Image is empty.')
else:
    plt.imshow(IMG)
    plt.show()
