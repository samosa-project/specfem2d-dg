
###############################################################
# HOW TO SEND AN EXAMPLE TO BATCH ON CALMIP.                  #
###############################################################

1) Make sure SPECFEM-DG is compiled with intel modules loaded. From the root folder, run:
    module purge; module load intel intelmpi; make clean; make all

2) From the EXAMPLES folder, run:
    ./prepare_example.sh EXAMPLENAME
  where EXAMPLENAME is the example's name.
  
  This script prepares the mesh and sets up symbolic links to the executables and to the parameter files.

3) From the EXAMPLES folder, run:
    ./send_example_to_execution_CALMIP.sh NNODES NCORES EXAMPLENAME TIMELIMIT
  where
    NNODES is the number of nodes,
    NCORES is the number of cores (NNODES * NCORES must match the number of processes specified in the parfile),
    EXAMPLENAME is the example's name, and
    TIMELIMIT is the time limit (format hours:minutes:seconds).
  
  This script submits the job to the Slurm controller.

###############################################################

