#include "PETScWrapper.h"

using namespace std;

//==============================================================================

int initPETScMat(Mat *A, int squareSize, int nonZeros)
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
  
}

//==============================================================================
