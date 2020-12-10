#ifndef PETSCWRAPPER_H
#define PETSCWRAPPER_H 

#include <petsc.h>
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std;

//==============================================================================

int initPETScMat(Mat *A, int squareSize, int nonZeros);
int initPETScRectMat(Mat *A, int rows, int cols, int nonZeros);
int initPETScVec(Vec *A, int size);
int eigenVecToPETScVec(Eigen::VectorXd *x_e,Vec *x_p);

//==============================================================================

#endif
