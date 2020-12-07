#ifndef PETSCWRAPPER_H
#define PETSCWRAPPER_H 

#include <petsc.h>

using namespace std;

//==============================================================================

void initPETScMat(Mat *A, double squareSize, double nonZeros);
void initPETScVec(Vec *A, double squareSize, double nonZeros);

//==============================================================================

#endif
