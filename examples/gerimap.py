#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def register():
    gm = [(0, 0, 0), (.15, .15, .5), (.3, .15, .75),
          (.6, .2, .50), (1, .25, .15), (.9, .5, 0),
          (.9, .75, .1), (.9, .9, .5), (1, 1, 1)]
    # Regular (dark-to-light)
    gerimap = LinearSegmentedColormap.from_list('GeriMap', gm)
    # Reversed (light-to-dark)
    gerimap_r = LinearSegmentedColormap.from_list('GeriMap_r', gm[::-1])

    plt.register_cmap(cmap=gerimap)
    plt.register_cmap(cmap=gerimap_r)

