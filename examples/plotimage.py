#!/usr/bin/env python3
#
# Loads a raw SOFT image and plots it.
# The first command-line argument gives the name of
# the SOFT output file to load. If a second argument
# is given, the generated figure is saved to a file
# with the name specified by the second argument.
# Otherwise, if no argument is given, the figure is
# shown on screen.
#
# ./plotimage.py infile [outfile]
#
# infile     -- Name of input file (i.e. SOFT output image)
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


def loadImage(filename):
    """
    Loads a pixel image from a SOFT image
    output file (in HDF5 or MATv7.3 format)

    filename: Name of file to load
    """
    image = None

    with h5py.File(filename, 'r') as f:
        image = f['image'][:,:].T

    return image


def plotImage(image):
    """
    Plot a SOFT image.

    image: Image matrix to plot.
    """
    global BACKGROUND_COLOR, COLORMAP

    fig = plt.figure()
    fig.patch.set_facecolor(BACKGROUND_COLOR)

    ax = fig.add_subplot(111)
    ax.imshow(image, cmap=COLORMAP, interpolation=None)
    ax.set_axis_off()

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

    fig.canvas.print_figure(outname, bbox_inches='tight', pad_inches=0, facecolor=fcolor, dpi=600)


def plothelp():
    print("./plotimage.py infile [outfile]")
    print("\nPlot a SOFT radiation image")
    print("infile     -- Name of input SOFT image (in HDF5/MATv7.3 format)")
    print("outfile    -- (optional) Name of output image. Filename extension")
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

    image   = loadImage(infile)
    ax, fig = plotImage(image)

    if outfile is not None:
        saveFigure(ax, fig, outfile)
    else:
        plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])

