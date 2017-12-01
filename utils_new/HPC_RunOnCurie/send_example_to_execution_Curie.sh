#!/bin/bash
# $1 = number of tasks
# $2 = example name
# $3 = time limit

if [[ $# -ne 3 ]]; then
echo
echo "> 2 arguments needed:"
echo "  - the number of tasks,"
echo "  - the example's name (a directory with that name must exist), and"
echo "  - the time limit (in seconds)."
echo "> Script will now exit."
echo
exit -1
fi

echo
echo "Tasks: $1."
echo "Example to be ran: '$2'."
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

cd ./$2/

# Stores output.
echo
echo "> Storing output (`date`)."
echo
cp ./DATA/*SOURCE* ./DATA/*STATIONS* ../../src/specfem2D/boundary_terms_DG.f90 ./OUTPUT_FILES

# Batch.
echo
echo "> Sending solver to queue (`date`)."
echo
ccc_msub -r $2 -n $1 -c 1 -T $3 -o %I.output -e %I.error -q standard -A gen10249 -x ../batch_example_job.sh
