#!/usr/bin/env python3
#
# Loads a raw SOFT Green's function of format '12'.
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


def loadGreen(filename):
    """
    Loads a Green's function output file
    (in HDF5 or MATv7.3 format)

    filename: Name of file to load
    """
    image = None

    tos = lambda v : "".join(map(chr, v[:,:][:,0].tolist()))

    with h5py.File(filename, 'r') as f:
        func = f['func'][:,:]
        par1 = f['param1'][:,:]
        par2 = f['param2'][:,:]

        frmt = tos(f['type'])
        par1n = tos(f['param1name'])
        par2n = tos(f['param2name'])

    # Verify format
    if frmt != '12':
        raise Exception("Invalid format of Green's function: '{0}'.".format(frm))

    # Reshape Green's function
    GF = np.reshape(func, [par1.size, par2.size]).T
    PPAR, PPERP = None, None

    if par1n == 'ppar' and par2n == 'pperp':
        PPAR, PPERP = par1, par2
    elif par1n == 'pperp' and par2n == 'ppar':
        PPAR, PPERP = par2, par1
        GF = GF.T
    else:
        raise Exception("Unrecognized momentum-space parameters in Green's function: '{0}' & '{1}'.".format(par1n, par2n))

    return PPAR, PPERP, GF


def plotGreen(ppar, pperp, gf):
    """
    Plot a SOFT Green's function.

    ppar:  Parallel momentum grid vector.
    pperp: Perpendicular momentum grid vector.
    gf:    Green's function.
    """
    global COLORMAP

    fig = plt.figure()

    PPAR, PPERP = np.meshgrid(ppar, pperp)

    ax = fig.add_subplot(111)
    F = gf * PPERP
    F = F / np.amax(F)
    ax.contourf(-PPAR, PPERP, F, cmap=COLORMAP)
    
    ax.set_xlabel('$p_\parallel / mc$')
    ax.set_ylabel('$p_\perp / mc$')

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
    print("./plotgreen12.py infile [outfile]")
    print("\nPlot a SOFT Green's function of format '12'")
    print("infile     -- Name of input SOFT Green's function (in HDF5/MATv7.3 format)")
    print("outfile    -- (optional) Name of output figure. Filename extension")
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

    PPAR, PPERP, GF = loadGreen(infile)
    ax, fig         = plotGreen(PPAR, PPERP, GF)

    if outfile is not None:
        saveFigure(ax, fig, outfile)
    else:
        plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])

