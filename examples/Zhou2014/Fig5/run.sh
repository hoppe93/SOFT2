#!/usr/bin/bash

source ../../config.sh 

###########
# Run SOFT
###########
mkdir -p data

time $SOFT a
time $SOFT b
time $SOFT c
time $SOFT d
time $SOFT e
time $SOFT f
time $SOFT g
time $SOFT h

###################
# Generate figures
###################
mkdir -p figures
mkdir -p topviews

# Orbits ((a) & (e))
../../plotorbits.py data/a-orbits.mat figures/a.pdf
../../plotorbits.py data/e-orbits.mat figures/e.pdf

# Images
../../plotimage.py data/b-image.mat figures/b.pdf
../../plotimage.py data/c-image.mat figures/c.pdf
../../plotimage.py data/d-image.mat figures/d.pdf

../../plotimage.py data/f-image.mat figures/f.pdf
../../plotimage.py data/g-image.mat figures/g.pdf
../../plotimage.py data/h-image.mat figures/h.pdf

# Top views
../../plottopview.py data/b-topview.mat topviews/b.pdf
../../plottopview.py data/c-topview.mat topviews/c.pdf
../../plottopview.py data/d-topview.mat topviews/d.pdf

../../plottopview.py data/f-topview.mat topviews/f.pdf
../../plottopview.py data/g-topview.mat topviews/g.pdf
../../plottopview.py data/h-topview.mat topviews/h.pdf

