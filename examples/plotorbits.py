#!/usr/bin/env python3
#
# Loads a set of SOFT orbits and plots them.
# The first command-line argument gives the name of
# the SOFT output file to load. If a second argument
# is given, the generated figure is saved to a file
# with the name specified by the second argument.
# Otherwise, if no argument is given, the figure is
# shown on screen.
#
# ./plottopview.py infile [outfile]
#
# infile     -- Name of input file (i.e. SOFT orbits file)
# outfile    -- Optional. If given, saves the generated
#               figure to a file with name 'outfile'.
#
# By: Mathias Hoppe, 2019
###############################

import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import gerimap
import matplotlib.ticker


# Configuration
FONTSIZE         = 16
LINEWIDTH        = 2
LINESTYLE        = 'k-'


def loadOrbits(filename):
    """
    Loads a set of orbits from a SOFT output
    file (in HDF5 or MATv7.3 format).

    filename: Name of file to load
    """
    T = None
    X = None

    with h5py.File(filename, 'r') as f:
        T = f['t'][:,:]
        X = f['x'][:,:]

    n, nt = T.shape

    X = np.reshape(X, (n, nt, 3))
    return T, X


def plotOrbits(X):
    """
    Plot a SOFT topview.

    X:    n-by-3-by-nT
          List of cartesian coordinates of each of
          the 'n' particles, in each of the 'nT'
          time steps.
    """
    global FONTSIZE, LINEWIDTH, LINESTYLE

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for x in X:
        R = np.sqrt(x[:,0]**2 + x[:,1]**2)
        Z = x[:,2]

        ax.plot(R, Z, LINESTYLE, linewidth=LINEWIDTH)

    ax.set_aspect('equal', 'box')
    ax.set_xlabel('Major radius R (m)', fontsize=FONTSIZE)
    ax.set_ylabel('Height Z (m)', fontsize=FONTSIZE)

    return ax, fig


def saveFigure(ax, fig, outname):
    """
    Save the plotted figure to a file.

    ax:      Matplotlib axis object.
    fig:     Matplotlib figure object.
    outname: Name of file to save figure to.
    """
    fig.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    fig.canvas.print_figure(outname, bbox_inches='tight', pad_inches=0)


def plothelp():
    print("./plottopview.py infile [outfile]")
    print("\nPlot a SOFT radiation topview.")
    print("infile     -- Name of input SOFT topview (in HDF5/MATv7.3 format)")
    print("outfile    -- (optional) Name of output topview. Filename extension")
    print("              determines the file type. Any file type supported by")
    print("              matplotlib can be used (common: PDF, PNG, EPS).")

def main(argv):
    infile  = None
    outfile = None

    if len(argv) >= 1:
        infile = argv[0]
    if len(argv) == 2:
        outfile = argv[1]
    if len(argv) > 2 or len(argv) == 0:
        print('ERROR: Too many command-line arguments given.')
        plothelp()
        sys.exit(-1)

    gerimap.register()

    T, X    = loadOrbits(infile)
    ax, fig = plotOrbits(X)

    if outfile is not None:
        saveFigure(ax, fig, outfile)
    else:
        plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])

