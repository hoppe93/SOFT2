#!/usr/bin/env python3
#
# Loads a raw SOFT top view and plots it.
# The first command-line argument gives the name of
# the SOFT output file to load. If a second argument
# is given, the generated figure is saved to a file
# with the name specified by the second argument.
# Otherwise, if no argument is given, the figure is
# shown on screen.
#
# ./plottopview.py infile [outfile]
#
# infile     -- Name of input file (i.e. SOFT topview output)
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
BACKGROUND_COLOR = 'black'
COLORMAP         = 'GeriMap'
LINECOLOR        = 'white'
LINEWIDTH        = 2


def loadTopview(filename):
    """
    Loads a topview from a SOFT topview
    output file (in HDF5 or MATv7.3 format)

    filename: Name of file to load
    """
    topview = None
    wall    = None

    with h5py.File(filename, 'r') as f:
        #detpos  = f['detectorPosition'][0,:].T
        detpos  = f['detectorPosition'][:].T
        topview = f['image'][:,:].T
        wall    = f['wall'][:,:]

    return topview, wall, detpos


def plotTopview(topview, wall, detpos):
    """
    Plot a SOFT topview.

    topview: Image matrix to plot.
    """
    global BACKGROUND_COLOR, COLORMAP, LINECOLOR, LINEWIDTH

    fig = plt.figure()
    fig.patch.set_facecolor(BACKGROUND_COLOR)

    # Get tokamak wall extents
    rwall = wall[0,:]
    rmin, rmax = np.amin(rwall), np.amax(rwall)
    extent = rmax*np.array([-1,1,-1,1])

    ax = fig.add_subplot(111)
    ax.imshow(topview, cmap=COLORMAP, interpolation=None, extent=extent)
    ax.set_axis_off()

    # Plot tokamak walls
    t = np.linspace(0, 2*np.pi, 100)
    sint, cost = np.sin(t), np.cos(t)

    xo, xi = rmax*cost, rmin*cost
    yo, yi = rmax*sint, rmin*sint

    ax.plot(xo, yo, linewidth=LINEWIDTH, color=LINECOLOR)
    ax.plot(xi, yi, linewidth=LINEWIDTH, color=LINECOLOR)

    # Plot detector position
    dp = detpos
    ax.plot(detpos[0], detpos[1], 'rx', markeredgewidth=5, markersize=14)

    return ax, fig


def saveFigure(ax, fig, outname):
    """
    Save the plotted figure to a file.

    ax:      Matplotlib axis object.
    fig:     Matplotlib figure object.
    outname: Name of file to save figure to.
    """
    ax.set_axis_off()
    fig.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)

    ax.get_xaxis().set_major_locator(matplotlib.ticker.NullLocator())
    ax.get_yaxis().set_major_locator(matplotlib.ticker.NullLocator())

    fcolor = fig.patch.get_facecolor()

    fig.canvas.print_figure(outname, bbox_inches='tight', pad_inches=0, facecolor=fcolor)


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

    topview, wall, detpos = loadTopview(infile)
    ax, fig               = plotTopview(topview, wall, detpos)

    if outfile is not None:
        saveFigure(ax, fig, outfile)
    else:
        plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])

