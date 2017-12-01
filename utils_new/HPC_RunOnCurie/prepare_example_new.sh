#!/bin/bash
# Runs mesher in serial.

parfileName="parfile_input"
sourcefileName="source_input"

if [[ $# -ne 1 ]]; then
echo
echo "> 1 argument needed:"
echo "  - the example's name (a directory with that name must exist)."
echo "> Script will now exit."
echo
exit -1
fi

simulationDir=$1 # This simulation's directory.

echo
echo "> Preparing example '$simulationDir' (`date`)."
echo
echo "> (This will take about 5 minutes.)"
echo

cd ./$simulationDir/
currentdir=`pwd`

# Sets up directory structure in current example directory.
echo
echo "> Setting up example's directory structure."
echo

mkdir -p OUTPUT_FILES
mkdir -p DATA

# Cleans current output and data files.
rm -rf DATA/*
rm -rf OUTPUT_FILES/*

# Sets up local DATA/ directory.
cd DATA/
ln -s ../$parfileName Par_file
ln -s ../$sourcefileName SOURCE
cd $currentdir

# Links executables.
rm -f xmeshfem2D xspecfem2D
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D

# Stores setup files.
cp DATA/Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# Save model configuration.
cp ../../src/specfem2D/compute_forces_acoustic_DG.f90 OUTPUT_FILES/compute_forces_acoustic_DG.f90
cp ../../src/specfem2D/boundary_terms_DG.f90 OUTPUT_FILES/boundary_terms_DG.f90

# Runs database generation.
# TODO: change this if an external mesh is to be used.
echo
echo "> Running mesher."
echo
./xmeshfem2D
