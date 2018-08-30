#!/bin/bash
# $1 = number of tasks
# $2 = example name
# $3 = time limit
# $4 = project number
# $5 = partition

if [[ $# -ne 5 ]]; then
  echo
  echo "> 5 arguments needed:"
  echo "  - the number of tasks,"
  echo "  - the example's name (a directory with that name must exist),"
  echo "  - the time limit (in seconds),"
  echo "  - the project number (##### in gen#####), and"
  echo "  - the partition (see ccc_mpinfo)."
  echo "> Script will now exit."
  echo
  exit -1
fi

echo
echo "> Example '$2' will be sent to execution, over $1 tasks, with $3 seconds time limit, on project 'gen$4'."
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
#echo "> [DEBUG] Command line: ccc_msub -r $2 -n $1 -c 1 -T $3 -o %I.output -e %I.error -q standard -A gen$4 -x ../batch_example_job.sh"
echo
ccc_msub -r $2 -n $1 -c 1 -T $3 -o %I.output -e %I.error -q $5 -A gen$4 -x ../batch_example_job.sh

# Check exit status.
exit_status=$?
if [ $exit_status -ne 0 ]; then
  echo "> [ERROR] Job submission failed."
  exit $exit_status
else
  echo "> Use 'ccc_mpp -u USERNAME' to check status of running jobs."
  echo "> Use 'ccc_mdel ID' to terminate job (where ID is obtained by checking status of running jobs)."
fi
