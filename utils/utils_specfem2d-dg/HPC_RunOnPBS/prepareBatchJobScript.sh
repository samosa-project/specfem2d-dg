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

JOBTAG="SPECFEM_JOB"
runThisExampleInFile="runThisExample.in"

curDir=$(pwd -P) # The "-P" option is to avoid symbolic links shenanigans.

if [[ ! -d ${curDir}/${example} ]]; then
echo "> [ERROR] EXAMPLE '$example' does not exist!"
exit -1
fi

if [[ ! -f ${curDir}/${runThisExampleInFile} ]]; then
echo "> ${runThisExampleInFile} not found."
exit -1
fi

# Move into EXAMPLE of interest.
cd "$curDir/$example/"

# Cleanup last run, if any.
if [[ -d "OUTPUT_FILES" ]]; then
id=$(ls ${JOBTAG}* | grep -oe "[0-9]\+" | head -n 1)
mv $slurmFile OUTPUT_FILES
mv ${JOBTAG}* OUTPUT_FILES
mv OUTPUT_FILES OUTPUT_FILES_$id
fi

# Write header.
echo "#!/bin/bash" > $slurmFile # THIS OVERWRITES ANY PREVIOUSLY EXISTING FILE
echo "#PBS -N ${JOBTAG}" >> $slurmFile # Appends.
echo "#PBS -q array-sn" >> $slurmFile # Appends. For some reason the array-sn queue works while the mpi-sn queue does not.
echo "#PBS -l select=$nnodes:ncpus=$ncores:mpiprocs=$ncores:mem=$ram" >> $slurmFile # Appends.
echo "#PBS -l walltime=$timelimit" >> $slurmFile # Appends.
echo "#PBS -m abe" >> $slurmFile # Appends.
echo "#PBS -M leo.martire@jpl.nasa.gov" >> $slurmFile # Appends.

# Load modules of interest.
echo "" >> $slurmFile # Appends.
echo "module load gcc/8.2.0 openmpi/gcc/64/4.0.4" >> $slurmFile # Appends.

# Add a command to explicitely cd into the EXAMPLE (since it seems PBS starts from wherever user lands when connecting).
echo "" >> $slurmFile # Appends.
echo "cd $curDir/$example/" >> $slurmFile # Appends.

# Copy-paste a custom "runThisExample.in" subscript in.
echo "" >> $slurmFile # Appends.
cat "$curDir/${runThisExampleInFile}" >> $slurmFile # Appends.

# Submit.
chmod u+x $slurmFile
qsub "$curDir/$example/$slurmFile"

# Move back.
cd $curDir

# Print queue.
ssh gattaca-edge-hn1 qstat -t -a -u lmartire 2>/dev/null; echo "To kill a job: qdel JOBID";echo "To show all jobs: qstat -a"
