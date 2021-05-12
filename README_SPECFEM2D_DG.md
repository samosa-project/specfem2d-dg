# SPECFEM2D-DG Configuration and Installation

## First Installation
### 1. Select the right Makefile for SCOTCH
  1. Check your architecture:
    ```
    uname -a
    ```
  2. Go into SCOTCH's src:
    ```
    cd ./src/meshfem2D/scotch/src
    ```
  3. Copy the appropriate Makefile from the Make.inc folder. E.g.:
    ```
    cp Make.inc/Makefile.inc.i686_mac_darwin10 Makefile.inc
    ```
  4. Go back to root:
    ```
    cd ../../../../
    ```

### 2. Configure SCOTCH
  1. Use the configurator:
    ```
    ./config_scotch.sh
    ```

  - Known error:
    ```
    library_random_f.c:75:3: error: implicit declaration of function 'SCOTCH_randomProc' is invalid in C99
          [-Werror,-Wimplicit-function-declaration]
      SCOTCH_randomProc (*procnum);
      ^
    ```
  - Solution:
    Edit `./src/meshfem2D/scotch/src/libscotch/scotch.h`.
    Add the following line to it (preferably around line 270, where it would make sense):
      `void SCOTCH_randomProc (int);`.
    Go into SCOTCH's src:
      `cd ./src/meshfem2D/scotch/src`.
    Re-make SCOTCH with this edit (do no re-run "./config_scotch.sh"):
      `make`.
    Go back to root and continue installation:
      `cd ../../../../`.
  - Cause:
    For some reason, ./config_scotch.sh replaces the content of ./src/meshfem2D/scotch/src/libscotch/scotch.h at each run, but with a file that does not have the explicit declaration for SCOTCH_randomProc. I could not find the source file used for that. If the source is found, simply add the line above to it and everything should work well.

  - Known error:
    ```
    undefined reference to `clock_gettime'
    ```
  - Solution:
    Check the glibc version installed on the system:
      `ldd --version`.
    If the version is < 2.17, a link to the real-time library might be needed. Add -lrt to the list of libraries in SCOTCH's Makefile. To do so, edit `./src/meshfem2D/scotch/src/Makefile.inc` and add "-lrt" at the end of the libraries' line ("LDFLAGS").
  - Cause:
    ??

  - Known error:
    ```
    /usr/bin/ld: cannot find -lz
    ```
  - Solution:
    Install the "zlib-dev" package.

  - Other tips:
    Check that the symbolic link "scotch" aims at the (last version) provided SCOTCH folder:
      `readlink -f scotch`.
    If not, make it be:
      `ln -s -T ./scotch_X.X.X ./scotch`

### 3. Configure SPECFEM
  1. Install an OpenMPI development package. Load its module. Check the module is indeed loaded ("module list" on Unix).
  2. Configure SPECFEM with GNU compilers by running:
    ```
    ./configure FC=gfortran CC=gcc MPIFC=mpif90 --with-mpi --enable-double-precision
    ```
  or with Intel compilers by running:
    ```
    ./configure FC=ifort CC=icc MPIFC=mpif90 --with-mpi --enable-double-precision
    ```
  
  - Known error:
    ```
    configure: error: MPI header not found; try setting MPI_INC.
    ```
  - Solution:
    Make sure you have an mpi package installed and loaded. Make sure you load the module corresponding to the compile you're using, e.g. with:
      `load openmpi/intel`.
  
  - Known error:
    ```
    configure: error: Fortran compiler cannot create executables
    ```
  - Solution:
    Make sure you have the right compilers installed (gcc, icc, gfortran, fort, etc.). On clusters, load the relevant package to have them available, e.g. with:
      `load intel/compiler`

### 4. Compile SPECFEM-DG
  1. Compile SPECFEM-DG:
    ```
    make clean; make all
    ```
  
  - Known error:
    ```
    undefined reference to `clock_gettime'"
    ```
  - Solution:
    Check the glibc version installed on the system:
      `ldd --version`.
    If the version is < 2.17, a link to the real-time library might be needed. Add -lrt to the list of libraries in SCOTCH's Makefile. To do so, edit `specfem-dg/src/meshfem2D/scotch/src/Makefile.inc` and add "-lrt" at the end of the libraries' line ("LDFLAGS").

  - Known error:
    ```
    Error: Array 'some_array' at (1) is larger than limit set by '-fmax-stack-var-size=', moved from stack to static storage. This makes the procedure unsafe when called recursively, or concurrently from multiple threads. Consider using '-frecursive', or increase the '-fmax-stack-var-size=' limit, or change the code to use an ALLOCATABLE array. [-Werror=surprising]
    ```
  Solution:
    Edit `./Makefile` (SPECFEM's main makefile).
    Add "-frecursive" to the list "FLAGS_CHECK".
    Recompile (do not re-run "./configure FC=..."):
      `make clean; make all`.

## Next Installations
  Now, SPECFEM2D-DG is configured for your machine. If you modify the source code, you only need to recompile it with:
  ```
    make
  ```
  or
  ```
    make clean; make all
  ```

# Running SPECFEM2D-DG Simulations
The simulations' method of execution varies whether you are on a local machine or on a cluster. On a local machine, running the script provided in the folders which are under "./EXAMPLES" should be enough. On a cluster, please refer to its documentation.

## Preparing a Run
Configure your own simulation case in the "./EXAMPLES/" folder.
1. Configure global parameters in the file "./EXAMPLES/my_simulation/PARFILE".
2. Adapt the mesh in the file "./EXAMPLES/my_simulation/INTERFACESFILE"
3. Configure source parameters in the file "./EXAMPLES/my_simulation/SOURCEFILE"
4. Note: the names of the files "PARFILE", "INTERFACESFILE", and "SOURCEFILE" may vary, depending of what's inside the "./EXAMPLES/my_simulation/run_this_example.sh" script.
5. Run the simulation using the "./EXAMPLES/my_simulation/run_this_example.sh" script.
