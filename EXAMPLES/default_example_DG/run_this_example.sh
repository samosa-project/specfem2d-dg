#!/bin/bash
# This script runs mesher and solver using the example setup.

parfile_input_name="parfile_input"
source_input_name="source_input"
interfacefileName="interfaces_input"

# Read the number of processes to use from the value of the "NPROC" line in the parfile.
NPROC=$(grep -e "NPROC *= *[0-9]*" $parfile_input_name | grep -oe "[0-9]*")

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
ln -s ../$parfile_input_name Par_file
ln -s ../$source_input_name SOURCE
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
cp ./DATA/Par_file ./OUTPUT_FILES/input_parfile
cp ./DATA/SOURCE ./OUTPUT_FILES/input_source
cp ./$interfacefileName OUTPUT_FILES/input_interfaces
cp ./atmospheric_model.dat ./OUTPUT_FILES/input_atmospheric_model.dat
cp ./external_bottom_forcing.dat ./OUTPUT_FILES/input_EBF.dat
# Save model config
cp ../../src/specfem2D/compute_forces_acoustic_DG.f90 OUTPUT_FILES/compute_forces_acoustic_DG.f90
cp ../../src/specfem2D/boundary_terms_DG.f90 OUTPUT_FILES/boundary_terms_DG.f90

# Runs database generation.
echo
echo ">> Starting mesher."
./xmeshfem2D

# Check if mesher ran well.
exit_status=$?
if [ $exit_status -ne 0 ] || [ ! -f ./OUTPUT_FILES/Database00000 ]; then
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
