#ifndef MultiPhysicsCoupledQD_H
#define MultiPhysicsCoupledQD_H

#include "Mesh.h"
#include "Materials.h"
#include "GreyGroupSolver.h"
#include "SingleGroupDNP.h"
#include "WriteData.h"

using namespace std;

class HeatTransfer;
class MultiGroupDNP;
class GreyGroupQD;

//==============================================================================
//! Contains precursor, heat, and grey group quasidiffusion objects.
//!   Also owns and is responsible for solving the coupled system each of
//!  of those object contribute to.

class MultiPhysicsCoupledQD
{
  public:
    // Constructor
    MultiPhysicsCoupledQD(Materials * myMats,\
        Mesh * myMesh,\
        YAML::Node * myInput);

    // Variables
    Eigen::SparseMatrix<double,Eigen::RowMajor> A;
    Eigen::VectorXd x,xPast,b;
    string outputDir = "MPQD/";
    double epsMPQD = 1E-6;
    int nUnknowns;
   
    // Functions 
    void fluxSource(int iZ,int iR,int iEq,double coeff,\
      Eigen::SparseMatrix<double,Eigen::RowMajor> * myA);
    void fluxSource(int iZ,int iR,int iEq,double coeff,\
      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> * myA);
    int dnpSource(int iZ,int iR,int iEq,double coeff,\
      Eigen::SparseMatrix<double,Eigen::RowMajor> * myA);
    void initializeXPast();
    void buildLinearSystem();
    void buildSteadyStateLinearSystem();
    void solveLinearSystem();
    void solveLinearSystemIterative(Eigen::VectorXd xGuess);
    int solveSuperLU();
    int solveIterativeDiag(Eigen::VectorXd xGuess);
    int solveIterativeILU(Eigen::VectorXd xGuess);
    void solveTransient();
    void updateVarsAfterConvergence();
    void updateSteadyStateVarsAfterConvergence();
    void writeVars();
    void printVars();
    void checkOptionalParams();
    int preconditioner = 1;

    // PETSc variables
    Vec x_p,xPast_p,b_p;
    Vec xPast_p_seq;
    Mat A_p;
    KSP ksp;
    PC pc;

    // Pointers
    HeatTransfer * heat;
    MultiGroupDNP * mgdnp;
    GreyGroupQD * ggqd;

  private:
    Materials * mats;
    Mesh * mesh;
    YAML::Node * input;
    const int iluPreconditioner = 0, diagPreconditioner = 1;

};

//==============================================================================

#endif
