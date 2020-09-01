###############################
# SPECFEM-DG CONFIGURATION    #
###############################

# First Installation. #########
1) Configure SCOTCH:
  Try:
    ./config_scotch.sh
  If it fails, we need to setup a bit more. Move to "./src/meshfem2D".
  
  Check that the symbolic link "scotch" aims at the (last version) provided SCOTCH folder:
    readlink -f scotch
  If not, make it be:
    ln -s -T ./scotch_X.X.X ./scotch

  Get back to "./specfem-dg". Try again:
    ./config_scotch.sh
  If it fails, we need to setup a bit more. Move to "./src/meshfem2D/scotch/src" in order to choose the appropriate Makefile. First, check your architecture:
    arch
  Copy the appropriate makefile from the "./Make.inc" folder. For example, on EOS (CALMIP) machines:
    cp ./Make.inc/Makefile.inc.x86-64_pc_linux2 Makefile.inc

  Error "undefined reference to `clock_gettime'": check the glibc version installed on the system with:
    ldd --version
  If the version is < 2.17, a link to the real-time library might be needed. Add -lrt to the list of libraries in SCOTCH's Makefile. To do so, edit "./src/meshfem2D/scotch/src/Makefile.inc" and add "-lrt" at the end of the libraries' line ("LDFLAGS").
  
  Error "/usr/bin/ld: cannot find -lz": install the "zlib-dev" package.
  
  If the scripts fails again, contact the SPECFEM-DG support.

2) Configure SPECFEM:
  Install an OpenMPI development package. Load its module. Check the module is indeed loaded:
    module list
  Configure SPECFEM with GNU compilers by running:
    ./configure FC=gfortran CC=gcc MPIFC=mpif90 --with-mpi --enable-double-precision
  or with Intel compilers by running:
    ./configure FC=ifort CC=icc MPIFC=mpif90 --with-mpi --enable-double-precision

3) Compile SPECFEM-DG:
  Compile SPECFEM-DG:
    make clean; make all
  Error "undefined reference to `clock_gettime'": check the glibc version installed on the system with:
    ldd --version
  If the version is < 2.17, a link to the real-time library might be needed. Add -lrt to the list of libraries in SCOTCH's Makefile. To do so, edit "specfem-dg/src/meshfem2D/scotch/src/Makefile.inc" and add "-lrt" at the end of the libraries' line ("LDFLAGS").

4) The simulations' method of execution varies whether you are on a local machine or on a cluster. On a local machine, running the script provided in the folders which are under "./EXAMPLE" should be enough. On a cluster, please report to its documentation.
###############################

# Running. ####################
1) Configure your own simulation: it happens in the "./EXAMPLES/" folder.
1.1) Configure global parameters in the file "./EXAMPLES/my_simulation/PARFILE".
1.2) Adapt the mesh in the file "./EXAMPLES/my_simulation/INTERFACESFILE"
1.3) Configure source parameters in the file "./EXAMPLES/my_simulation/SOURCEFILE"
1.4) Note: the names of the files "PARFILE", "INTERFACESFILE", and "SOURCEFILE" may vary, depending of what's inside the "./EXAMPLES/my_simulation/run_this_example.sh" script.
1.5) Run the simulation using the "./EXAMPLES/my_simulation/run_this_example.sh" script.

2) [DEVELOPEMENT] How to configure DG FNS (full Navier-Stokes) boundary conditions.
2.1) Modify the source file "./src/specfem2d/boundary_terms.f90" according to the wanted boundary conditions.
2.2) Recompile SPECFEM-DG:
       make
     or:
       make clean; make all
###############################

# Runtime errors. #############
*) Issue:    ###
   Cause:    ###
   Solution: ###
###############################