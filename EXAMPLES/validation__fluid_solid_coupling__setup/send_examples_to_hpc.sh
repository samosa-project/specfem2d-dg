#!/bin/bash

# EOS no VPN
remote=martire@olympe.calmip.univ-toulouse.fr


dest="${remote}:/tmpdir/martire/specfem-dg-master/EXAMPLES/"
source="/home/l.martire/Documents/SPECFEM/specfem-dg-master/EXAMPLES/"
prefix="validation__fluid_solid_coupling__*"

echo " "
printf %"$(tput cols)"s |tr " " "*" # pretty line to separate script from other commands
printf %"$(tput cols)"s |tr " " "*"
curDir=$(pwd) # save current dir

echo "  Considering EXAMPLE files in ${source}."
echo "  Sending to remote at ${dest}."
echo " "

echo "  Building list of files to send."
cd $source
ls ${prefix}/parfile_input ${prefix}/interfaces_input ${prefix}/source_input > $curDir/filesToSend

echo "  Files to send:"
cat $curDir/filesToSend

echo "  Performing rsync."
rsync -v -e "ssh " -avz --progress ${source} --files-from=$curDir/filesToSend ${dest}

printf %"$(tput cols)"s |tr " " "*" # pretty line to separate script from other commands
printf %"$(tput cols)"s |tr " " "*"
echo " "

cd $curDir # return to current dir
