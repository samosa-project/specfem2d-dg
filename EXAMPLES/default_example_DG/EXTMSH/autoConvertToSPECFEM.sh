#!/bin/bash
# Author:        Léo Martire.
# Mail:          leo.martire@outlook.com
# Description:   GMSH meshing pipeline for SPECFEM2D.
# Last modified: See file metadata.
# Usage:         1) If called with 'auto' as argument, will automatically mesh and convert to SPECFEM standards.
#                2) If called with another or no argument, will open GMSH GUI enabling the user to mesh as they want, then convert to SPECFEM standards on closing of the GUI.
# Notes:         N/A.

# Save the current directory and script directory.
callingDir=$(pwd)
thisScriptDir=$(dirname $(realpath $0))

###############################################
# Converting .geo to .msh                     #
###############################################
cd $thisScriptDir
if [[ "$1" == auto ]]; then
  echo "[$0] Automatic meshing."
  gmsh $thisScriptDir/*.geo -2
else
  echo "[$0] Manual meshing."
  gmsh -open $thisScriptDir/*.geo
fi
  

###############################################
# Eventually open GMSH to show mesh.          #
###############################################
if [[ "$1" == auto ]]; then
  meshToShow="extMesh"
  echo "[$0] GMSH opened to show $meshToShow."
  gmsh -merge $thisScriptDir/$meshToShow.geo $thisScriptDir/$meshToShow.msh&
fi

###############################################
# Converting mesh to SPECFEM standards.       #
###############################################
echo "[$0] Locating GMSH->SPECFEM script."
# Convert script.
scriptName="LibGmsh2Specfem_convert_Gmsh_to_Specfem2D_official.py"
# Go through the folders containing this folder to find the conversion script.
found=-1
cd $thisScriptDir
while [ $found -ne 0 ]
do
  cd ..
  fullPathToScript=$(tree -fiP "$scriptName" | grep "$scriptName")
  found=$?
done
# Save the path to the conversion script.
absolutePathToScript="$(pwd)/$fullPathToScript"
echo "[$0] Converting to SPECFEM standards."
# Ask user for boundaries' types input.
echo "> Boundaries?"
echo "  > Top?"
read bt
echo "  > Left?"
read bl
echo "  > Bottom?"
read bb
echo "  > Right?"
read br
# Convert user's inputs to uppercase.
bt=$(echo $bt | tr a-z A-Z)
bl=$(echo $bl | tr a-z A-Z)
bb=$(echo $bb | tr a-z A-Z)
br=$(echo $br | tr a-z A-Z)
# Go back to the initial directory.
cd $thisScriptDir
# Run the conversion script.
python $absolutePathToScript *.msh -t $bt -l $bl -b $bb -r $br

# Eventually go back to calling directory for other instructions.
cd $callingDir
