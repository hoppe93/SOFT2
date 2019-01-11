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

###################
# Generate figures
###################
mkdir -p figures
mkdir -p topviews

# Orbits ((a) & (e))
../../plotorbits.py data/a-orbits.mat figures/a.pdf
../../plotorbits.py data/c-orbits.mat figures/c.pdf

# Images
../../plotimage.py data/b-image.mat figures/b.pdf
../../plotimage.py data/d-image.mat figures/d.pdf

# Top views
../../plottopview.py data/b-topview.mat topviews/b.pdf
../../plottopview.py data/d-topview.mat topviews/d.pdf

