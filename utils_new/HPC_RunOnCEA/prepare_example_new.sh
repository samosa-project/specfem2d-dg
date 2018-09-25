#!/bin/bash
# Runs mesher in serial.

parfileName="parfile_input"
sourcefileName="source_input"
interfacefileName="interfaces_input"

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

# Sets up directory structure in current example directoy.
echo
echo ">> Make 'DATA' and 'OUTPUT_FILES' folder."
mkdir -p OUTPUT_FILES
mkdir -p DATA

# Cleans current output and data files.
echo
echo ">> Clean the 'OUTPUT_FILES' folder."
rm -rf DATA/*
rm -rf OUTPUT_FILES/*

# Sets up local DATA/ directory.
echo
echo ">> Make symbolic links to parfile and sourcefile in the 'DATA' folder."
cd DATA/
ln -s ../$parfileName Par_file
ln -s ../$sourcefileName SOURCE
cd $currentdir

# Links executables.
echo
echo ">> Make symbolic links to executables."
rm -f xmeshfem2D xspecfem2D
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D

# Stores input files.
echo
echo ">> Store setup to the 'OUTPUT_FILES' folder."
cp ./DATA/Par_file ./OUTPUT_FILES/input_parfile
cp ./DATA/SOURCE ./OUTPUT_FILES/input_source
cp ./$interfacefileName ./OUTPUT_FILES/input_interfaces
cp ./atmospheric_model.dat ./OUTPUT_FILES/input_atmospheric_model.dat
if [ -e "external_bottom_forcing.dat" ]; then
  tar -czf ./OUTPUT_FILES/input_EBF.tar ./external_bottom_forcing.dat
fi
if [ -d "EXTMSH" ]; then
  tar -czf ./OUTPUT_FILES/input_EXTMSH.tar ./EXTMSH/*
fi
# Save model configuration.
cp ../../src/specfem2D/compute_forces_acoustic_DG.f90 OUTPUT_FILES/compute_forces_acoustic_DG.f90
cp ../../src/specfem2D/boundary_terms_DG.f90 OUTPUT_FILES/boundary_terms_DG.f90

# Runs database generation.
# TODO: change this if an external mesh is to be used.
echo
echo "> Running mesher."
echo
./xmeshfem2D
