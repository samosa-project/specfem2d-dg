
###############################################################
# HOW TO SEND AN EXAMPLE TO BATCH ON PBS-RUNNING MACHINES.    #
###############################################################

1) From the EXAMPLES folder, run:
    ./prepareBatchJobScript.sh NNODES NTASKS EXAMPLENAME TIMELIMIT RAM
  where
    NNODES is the number of nodes,
    NTASKS is the number of tasks (NNODES * NCORES must match the number of processes specified in the parfile),
    EXAMPLENAME is the example's name,
    TIMELIMIT is the time limit (under format hours:minutes:seconds), and
    RAM is the RAM per node (under format 256mb, 8gb, etc.).

  This script prepares the mesh, sets up symbolic links to the executables and to the parameter files, and submits the job to the PBS controller.

###############################################################

