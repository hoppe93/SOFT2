#!/usr/bin/bash

source ../config.sh 

###########
# Run SOFT
###########
mkdir -p data

time $SOFT gc
time $SOFT particle

###################
# Generate figures
###################
../plotorbits.py data/guiding-center.mat data/guiding-center.pdf
../plotorbits.py data/particle.mat data/particle.pdf

