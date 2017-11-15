
###############################################################
# HOW TO SEND AN EXAMPLE TO BATCH ON curie.                   #
###############################################################

1) Run:
    ./prepare_example.sh EXAMPLENAME
  where EXAMPLENAME is the example's name.

2) Edit "batch_example.slurm" and set the parameters.

3) Run:

4) Run:
    current_example=CURRENT_EXAMPLE_NAME
    cp .$current_example/DATA/*SOURCE* .$current_example/DATA/*STATIONS* ../src/specfem2D/boundary_terms_DG.f90 .$current_example/OUTPUT_FILES
  to save source information, stations, and boundary terms source file.

###############################################################

