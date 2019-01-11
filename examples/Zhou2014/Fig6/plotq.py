#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys


# Number of points to evaluate
# safety factor in
N    = 100
# Maximum minor radius
RMAX = 45


def plotqprofile(qtype, qa1, qa2=0):
    """
    Plot a safety factor profile.

    a: Normalized radius.
    q: Safety factor.
    """
    global N, RMAX

    outname = None
    argv = sys.argv[1:]

    if len(argv) == 1:
        outname = argv[0]
    elif len(argv) > 1:
        print('ERROR: Too many command-line arguments.')
        print('Expected, at most, name of file to save figure to.')
        sys.exit(-1)

    a = np.linspace(0, 1, N)
    q = None

    if qtype == 'constant':
        q = qa1*np.ones(a.shape)
    elif qtype == 'linear':
        q = qa1*a + qa2
    elif qtype == 'quadratic':
        q = qa1*a**2 + qa2
    elif qtype == 'exponential':
        q = np.exp(qa1) + qa2
    else:
        print("ERROR: Unrecognized safety factor type: '{0}'".format(qtype))

    fig = plt.figure(figsize=(5,4))
    ax  = fig.add_subplot(111)

    ax.plot(a*RMAX, q, 'b', linewidth=2)
    ax.set_xlabel('r (cm)', fontsize=16)
    ax.set_ylabel('q', fontsize=16)

    ax.set_xlim([0, RMAX])
    ax.set_ylim([0, 5])

    if outname is not None:
        fig.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        fig.canvas.print_figure(outname, bbox_inches='tight', pad_inches=0)
    else:
        plt.show()

