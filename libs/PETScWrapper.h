#ifndef PETSCWRAPPER_H
#define PETSCWRAPPER_H 

#include <petsc.h>

using namespace std;

//==============================================================================

int initPETScMat(Mat *A, int squareSize, int nonZeros);
int initPETScVec(Vec *A, int size);

//==============================================================================

#endif
