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
   
    // Functions 
    void fluxSource(int iZ,int iR,int iEq,double coeff);
    void dnpSource(int iZ,int iR,int iEq,double coeff);
    void initializeXPast();
    void buildLinearSystem();
    void solveLinearSystem();
    void solveLinearSystemParallel();
    void solveTransient();
    void updateVarsAfterConvergence();
    void writeVars();
    void printVars();
    void checkOptionalParams();

    // Pointers
    HeatTransfer * heat;
    MultiGroupDNP * mgdnp;
    GreyGroupQD * ggqd;

  private:
    Materials * mats;
    Mesh * mesh;
    YAML::Node * input;

};

//==============================================================================

#endif
