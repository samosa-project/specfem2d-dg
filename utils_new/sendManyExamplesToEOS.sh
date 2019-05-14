#!/bin/bash
# Author:        Léo Martire.
# Description:   Sends many parameter files to EOS.
# Last modified: See file metadata.
# Usage:         N/A.
# Notes:         N/A.


localEXAMPLEDirPath='/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/'
remoteEXAMPLEDirPath='/tmpdir/martire/specfem-dg-master/EXAMPLES/'
listOfSubFilesToSend="parfile_input,interfaces_input,source_input,run_this_example.sh,atmospheric_model.dat"

cd ${localEXAMPLEDirPath}

nargneed=1
if [[ $# -ne $nargneed ]]; then
  echo
  echo "> This script needs $nargneed argument, a list of EXAMPLE folders separated by commas, e.g.:"
  echo "    default_example_DG,anotherOne,wHaTeVeR"
  echo "> Script will now exit."
  echo
  exit -1
fi

echo " "
echo " "
echo " "

listOfEXAMPLES=$1
#listOfEXAMPLES='test1,test2'

command="rsync -avz -e 'ssh' --relative {$listOfEXAMPLES}/{$listOfSubFilesToSend} martire@olympe.calmip.univ-toulouse.fr:${remoteEXAMPLEDirPath}"
echo "Running command '$command'."
echo " "

eval $command

echo " "
echo " "
echo " "

#rsync -avz --relative test1/{$listOfSubFilesToSend} martire@olympe.calmip.univ-toulouse.fr:${remoteEXAMPLEDirPath}

#rsync -avz -e 'ssh' '/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/LNSABC_info/README' martire@olympe.calmip.univ-toulouse.fr:${remoteEXAMPLEDirPath}
