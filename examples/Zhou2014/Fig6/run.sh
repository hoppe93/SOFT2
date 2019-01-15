#!/usr/bin/bash

source ../../config.sh

###########
# Run SOFT
###########
mkdir -p data

time $SOFT e
time $SOFT f
time $SOFT g
time $SOFT h

###################
# Generate figures
###################
mkdir -p figures
mkdir -p topviews

# Safety factor profiles
./gena.py figures/a.pdf
./genb.py figures/b.pdf
./genc.py figures/c.pdf
./gend.py figures/d.pdf

# Images
../../plotimage.py data/e-image.mat figures/e.pdf
../../plotimage.py data/f-image.mat figures/f.pdf
../../plotimage.py data/g-image.mat figures/g.pdf
../../plotimage.py data/h-image.mat figures/h.pdf

# Top views
../../plottopview.py data/e-topview.mat topviews/e.pdf
../../plottopview.py data/f-topview.mat topviews/f.pdf
../../plottopview.py data/g-topview.mat topviews/g.pdf
../../plottopview.py data/h-topview.mat topviews/h.pdf

