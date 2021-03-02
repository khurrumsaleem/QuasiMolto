#!/bin/bash

# Armadillo
cd armadillo/
cmake . -DCMAKE_INSTALL_PREFIX:PATH=.
make
make install
# SuperLU
cd ../superlu/
cmake . -DCMAKE_INSTALL_PREFIX:PATH=.
make
make install
# PETSc
cd ../petsc
./configure --with-fc=0
make PETSC_DIR=$(pwd) PETSC_ARCH=arch-linux-c-debug all
make PETSC_DIR=$(pwd) PETSC_ARCH=arch-linux-c-debug check
