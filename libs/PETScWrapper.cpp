#include "PETScWrapper.h"

using namespace std;

//==============================================================================

int initPETScMat(Mat * A, int squareSize, int nonZeros)
{
  PetscErrorCode ierr; 

  /* Initialize matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,A);
  CHKERRQ(ierr);
  
  /* Set size of matrix */
  ierr = MatSetSizes(*A,PETSC_DECIDE,PETSC_DECIDE,squareSize,squareSize);
  CHKERRQ(ierr);
  
  /* Pull in command line options */
  ierr = MatSetFromOptions(*A);
  CHKERRQ(ierr);
  
  /* Preallocate matrix memory */
  ierr = MatMPIAIJSetPreallocation(*A,nonZeros,NULL,nonZeros,NULL);
  CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(*A,nonZeros,NULL);
  CHKERRQ(ierr);

  return ierr;

}

int initPETScRectMat(Mat *A, int rows, int cols, int nonZeros)
{
  PetscErrorCode ierr; 

  /* Initialize matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,A);
  CHKERRQ(ierr);
  
  /* Set size of matrix */
  ierr = MatSetSizes(*A,PETSC_DECIDE,PETSC_DECIDE,rows,cols);
  CHKERRQ(ierr);
  
  /* Pull in command line options */
  ierr = MatSetFromOptions(*A);
  CHKERRQ(ierr);
  
  /* Preallocate matrix memory */
  ierr = MatMPIAIJSetPreallocation(*A,nonZeros,NULL,nonZeros,NULL);
  CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(*A,nonZeros,NULL);
  CHKERRQ(ierr);
  
  return ierr;

}


int initPETScVec(Vec *x, int size)
{
  PetscErrorCode ierr;
  
  /* Initialize matrix */
  ierr = VecCreate(PETSC_COMM_WORLD,x);
  CHKERRQ(ierr);
  
  /* Set size of vector */
  ierr = VecSetSizes(*x,PETSC_DECIDE,size);
  CHKERRQ(ierr);
  
  /* Pull in command line options */
  ierr = VecSetFromOptions(*x);
  CHKERRQ(ierr);
  
  return ierr;
  
}

int eigenVecToPETScVec(Eigen::VectorXd *x_e,Vec *x_p)
{
  PetscErrorCode ierr;
  
  for (int idx = 0; idx < (*x_e).size(); idx++)
  {
    ierr = VecSetValue(*x_p,idx,(*x_e)(idx),INSERT_VALUES);CHKERRQ(ierr); 
  }
  ierr = VecAssemblyBegin(*x_p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*x_p);CHKERRQ(ierr);
  
  return ierr;
}

int petscVecToEigenVec(Vec *x_p,Eigen::VectorXd *x_e)
{
  PetscErrorCode ierr;
  PetscScalar value;
  PetscInt vecSize;
  VecScatter     ctx;
  Vec temp;

  // Collect values from all procs
  VecScatterCreateToAll(*x_p,&ctx,&temp);
  VecScatterBegin(ctx,*x_p,temp,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,*x_p,temp,INSERT_VALUES,SCATTER_FORWARD);

  // Set set of Eigen vector
  ierr = VecGetSize(*x_p, &vecSize);
  (*x_e).setZero(vecSize);

  for (int idx = 0; idx < vecSize; idx++)
  {
    ierr = VecGetValues(temp,1,&idx,&value);CHKERRQ(ierr);
    (*x_e)(idx) = value;
  }
  
  VecScatterDestroy(&ctx);
  
  return ierr;
}


//==============================================================================
