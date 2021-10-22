# QuasiMolto [![Build Status](https://travis-ci.org/aaronjamesreynolds/QuasiMolto.svg?branch=master)](https://travis-ci.org/aaronjamesreynolds/QuasiMolto)

A multigroup, multiphysics, multilevel transport solver for circulating fuel reactor kinetics. 
Contains transient and steady-state solvers. 

## Build instructions

  * Set compilers (note: compatability confirmed for the 7.5.0 versions of `mpicc` and `mpicxx` distributed by GNU) 

    ```export CC=mpicc```
    
    ```export CXX=mpicxx```

  * Clone the QuasiMolto repository and its submodules

    ```git clone --recurse-submodules git@github.com:aaronjamesreynolds/QuasiMolto.git```

  * Build the third-party libraries
  
    ```cd QuasiMolto/TPLs/```
	
    ```chmod +x build-tpls.sh```
    
    ```source ./build-tpls.sh```
   
  * Make a build directory and compile `QuasiMolto`
  
    ```cd ../../../```
    
    ```mkdir build_directory```
    
    ```cd build_directory```
    
    ```cmake -DCMAKE_BUILD_TYPE=Release ../QuasiMolto/```
    
    ```make all```
    
    ```make test```
    
  * The `QuasiMolto` executable is in the `bin` directory within `build_directory`. 
    A simple input file is located in `QuasiMolto/examples/homogeneous/`. 
    A simulation directory can be built, and the sample input run, with the following commands. 
    
    ```cd ../```
    
    ```mkdir myRun```
    
    ```cd myRun```
    
    ```cp ../QuasiMolto/examples/homogeneous/* .```
    
    ```../build_directory/bin/QuasiMolto homogeneous.yaml```

![alt text](https://vignette.wikia.nocookie.net/monstermovies/images/4/46/Quasimodo.png/revision/latest?cb=20140628171627 "Quasi Moto")
