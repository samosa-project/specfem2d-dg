#!/bin/bash

slurmFile="batchJobScript.sh"

if [[ $# -ne 5 ]]; then
echo
echo "> 5 arguments needed:"
echo "  - the number of nodes,"
echo "  - the number of cores,"
echo "  - the example's name (a directory with that name must exist),"
echo "  - the time limit (under format hours:minutes:seconds), and"
echo "  - the RAM per node (under format 256mb, 8gb, etc.)."
echo "> Script will now exit."
echo
exit -1
fi
nnodes=$1
ncores=$2
example=$3
timelimit=$4
ram=$5

if [[ ! -d $example ]]; then
echo "> [ERROR] EXAMPLE '$example' does not exist!"
exit -1
fi

# Move into EXAMPLE of interest.
curDir=$(pwd)
cd "$curDir/$example/"

# Cleanup last run, if any.
if [[ -d "OUTPUT_FILES" ]]; then
id=$(ls PBSJOB* | grep -oe "[0-9]\+" | head -n 1)
mv $slurmFile OUTPUT_FILES
mv PBSJOB* OUTPUT_FILES
mv OUTPUT_FILES OUTPUT_FILES_$id
fi

# Write header.
echo "#!/bin/bash" > $slurmFile # THIS OVERWRITES ANY PREVIOUSLY EXISTING FILE
echo "#PBS -N PBSJOB" >> $slurmFile # Appends.
echo "#PBS -q mpi" >> $slurmFile # Appends.
echo "#PBS -l select=$nnodes:ncpus=$ncores:mpiprocs=$ncores:mem=$ram" >> $slurmFile # Appends.
echo "#PBS -l walltime=$timelimit" >> $slurmFile # Appends.

# Add a command to explicitely cd into the EXAMPLE (since it seems PBS starts from wherever user lands when connecting).
echo "" >> $slurmFile # Appends.
echo "cd $curDir/$example/" >> $slurmFile # Appends.

# Copy-paste a custom "runThisExample" subscript in.
echo "" >> $slurmFile # Appends.
cat "$curDir/runThisExample.in" >> $slurmFile # Appends.

# Submit.
chmod u+x $slurmFile
qsub "$curDir/$example/$slurmFile"

# Move back.
cd $curDir
