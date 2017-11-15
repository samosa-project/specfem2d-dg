#!/bin/bash
#MSUB -r MARTIRE_TEST_RUN # Job name.
#MSUB -n 4 # Number of tasks (product of the two following quantities).
#MSUB -c 2 # Number of cores.
#MSUB -N 2 # Number of nodes.
#MSUB -T 60 # Time limit (in seconds).
#MSUB -A gen10249 # Project ID.

example_name=full_D

set -x

cd ${BRIDGE_MSUB_PWD}

ccc_mprun xspecfem2D
