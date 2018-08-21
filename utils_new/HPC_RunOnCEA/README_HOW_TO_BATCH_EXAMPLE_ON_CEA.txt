
###############################################################
# HOW TO SEND AN EXAMPLE TO BATCH ON CEA MACHINES.            #
###############################################################

1) From the EXAMPLES folder, run:
    ./prepare_example.sh EXAMPLENAME
  where EXAMPLENAME is the example's name.

  This script prepares the mesh and sets up symbolic links to the executables and to the parameter files.

2) From the EXAMPLES folder, run:
    ./send_example_to_execution_CEA.sh NTASKS EXAMPLENAME TIMELIMIT PROJECTNUMBER PARTITION
  where
    NTASKS is the number of tasks (must correspond to NPROC in parfile),
    EXAMPLENAME is the example's name,
    TIMELIMIT is the time limit (in seconds),
    PROJECTNUMBER is the project number (##### in gen#####), and
    PARTITION is the partition (see ccc_mpinfo).
  
  This script submits the job to the Slurm controller.

###############################################################

