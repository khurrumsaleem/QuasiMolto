#ifndef PETSCWRAPPER_H
#define PETSCWRAPPER_H 

#include <petsc.h>
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"
#include <iostream>

using namespace std;

//==============================================================================

int initPETScMat(Mat * A, int squareSize, int nonZeros);
int initPETScRectMat(Mat * A, int rows, int cols, int nonZeros);
int initPETScVec(Vec * A, int size);
int eigenVecToPETScVec(Eigen::VectorXd * x_e,Vec * x_p);
int petscVecToEigenVec(Vec * x_p,Eigen::VectorXd * x_e);

//==============================================================================

#endif
