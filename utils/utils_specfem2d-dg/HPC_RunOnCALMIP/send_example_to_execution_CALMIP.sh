#!/bin/bash
# $1 = nombre de noeuds
# $2 = nombre de coeurs
# $3 = nom de l'exemple
# Quantity $1 x $2 corresponds to the wanted number of CPUs.

slurmFile="batch_example.slurm"

if [[ $# -ne 4 ]]; then
echo
echo "> 4 arguments needed:"
echo "  - the number of nodes,"
echo "  - the number of cores,"
echo "  - the example's name (a directory with that name must exist), and"
echo "  - the time limit (format hours:minutes:seconds)."
echo "> Script will now exit."
echo
exit -1
fi

result_nodes=$(( $1*$2 ))
echo
echo "Tasks: $result_nodes ($1 nodes * $2 cores)."
echo "Example to be ran: '$3'."
echo

# Loads the modules necessary to the execution of the script (and only those).
# echo
# echo "> Loading modules."
# echo
# module purge
# #module load intel/14.0.2.144 intelmpi/4.1.3.049
# module load openmpi
# module list
# export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

cd ./$3/

# Stores output.
echo
echo "> Storing output (`date`)."
echo
cp ./DATA/*SOURCE* ./DATA/*STATIONS* ../../src/specfem2D/boundary_terms_DG.f90 ./OUTPUT_FILES

# Batch.
echo
echo "> Sending solver to queue (`date`)."
echo
sbatch --nodes=$1 --ntasks-per-node=$2 --ntasks=$result_nodes --job-name=$3 --time=$4 ../$slurmFile
