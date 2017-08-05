#!/bin/bash
#

echo "Condiguration of SCOTCH"
currentdir=`pwd`

echo
echo "(will take about 5 minutes)"
echo

# sets up link to scotch
cd src/meshfem2D/
ln -s scotch_6.0.4 scotch

# Compile scotch
cd scotch/src/
make clean
make

echo
echo "done"
date
echo
