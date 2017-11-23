#!/bin/bash

# Loads the modules necessary to the execution of the script (and only those).
module purge
#module load intel/14.0.2.144 intelmpi/4.1.3.049
module load intel intelmpi
#module load openmpi
module list
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# Runs simulation.
srun xspecfem2D
#mpirun xspecfem2D
