#!/bin/bash

scriptName="LibGmsh2Specfem_convert_Gmsh_to_Specfem2D_official.py"

# Save the current directory.
thisScriptDir=$(dirname $(realpath $0))

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
