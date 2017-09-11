#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 5 minutes)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p OUTPUT_FILES
mkdir -p DATA

# sets up local DATA/ directory
cd DATA/
ln -s ../Par_file_fluid_solid Par_file
ln -s ../SOURCE_fluid_solid SOURCE
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

cd $currentdir

# links executables
rm -f xmeshfem2D xspecfem2D
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# Save model config
cp ../../src/specfem2D/compute_forces_acoustic_DG.f90 OUTPUT_FILES/compute_forces_acoustic_DG.f90
cp ../../src/specfem2D/boundary_terms.f90 OUTPUT_FILES/boundary_terms.f90

# runs database generation
echo
echo "  running mesher..."
echo
./xmeshfem2D
mpirun -np 4 ./xspecfem2D

# runs simulation
echo
echo "  running solver..."
echo
#./xspecfem2D

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
