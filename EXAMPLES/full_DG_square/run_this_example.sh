#!/bin/bash
# This script runs mesher and solver (in serial) using this example setup.

# Set this value to the value of the corresponding line in the parfile.
NPROC=4

echo
echo ">> Running example (`date`)."
currentdir=`pwd`

# Sets up directory structure in current example directoy.
echo
echo ">> Make 'DATA' and 'OUTPUT_FILES' folder."
mkdir -p OUTPUT_FILES
mkdir -p DATA

# Sets up local "DATA/" directory.
echo
echo ">> Make symbolic links to parfile and sourcefile in the 'DATA' folder."
cd DATA/
rm -f Par_file SOURCE
ln -s ../parfile_input Par_file
ln -s ../source_input SOURCE
cd ../

# Cleans output files.
echo
echo ">> Clean the 'OUTPUT_FILES' folder."
rm -rf OUTPUT_FILES/*

cd $currentdir

# Links executables.
echo
echo ">> Make symbolic links to executables."
rm -f xmeshfem2D xspecfem2D
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D

# Stores setup.
echo
echo ">> Store setup to the 'OUTPUT_FILES' folder."
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/
# Save model config
cp ../../src/specfem2D/compute_forces_acoustic_DG.f90 OUTPUT_FILES/compute_forces_acoustic_DG.f90
cp ../../src/specfem2D/boundary_terms_DG.f90 OUTPUT_FILES/boundary_terms_DG.f90

# Runs database generation.
echo
echo ">> Starting mesher."
./xmeshfem2D

# Check if mesher ran well.
exit_status=$?
if [ $exit_status -ne 0 ]; then
  echo
  echo ">> [ERROR] Mesher did not finish well. Exiting."
  exit $exit_status
fi

# Save source and stations info.
echo
echo ">> Finish storing output to the 'OUTPUT_FILES' folder."
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

# Runs simulation.
echo
echo ">> Starting solver."
mpirun -np $NPROC ./xspecfem2D

# Check if solver ran well.
exit_status=$?
if [ $exit_status -ne 0 ]; then
  echo
  echo ">> [ERROR] Solver did not finish well. Exiting."
  exit $exit_status
fi

echo
echo ">> Done (`date`). See results in the 'OUTPUT_FILES' folder."
