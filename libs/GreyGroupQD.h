#ifndef GreyGroupQD_H
#define GreyGroupQD_H

#include "Mesh.h"
#include "Materials.h"
#include "MultiGroupQD.h"
#include "SingleGroupQD.h"

using namespace std;

class MultiPhysicsCoupledQD;
class GreyGroupSolverBase;
class GreyGroupSolver;

//==============================================================================
//! Contains information and builds linear system for a grey group qd

class GreyGroupQD
{
  public:

    // VARIABLES

    int indexOffset,nUnknowns,nCurrentUnknowns;

    Eigen::MatrixXd sFlux;
    Eigen::MatrixXd sFluxR;
    Eigen::MatrixXd sFluxZ;
    Eigen::MatrixXd currentR;
    Eigen::MatrixXd currentZ;
    Eigen::MatrixXd sFluxPrev;
    Eigen::MatrixXd q;

    // Eddington factors
    Eigen::MatrixXd Err,ErrPrev;
    Eigen::MatrixXd Ezz,EzzPrev;
    Eigen::MatrixXd Erz,ErzPrev;
    Eigen::MatrixXd ErrAxial;
    Eigen::MatrixXd EzzAxial;
    Eigen::MatrixXd ErzAxial;
    Eigen::MatrixXd ErrRadial;
    Eigen::MatrixXd EzzRadial;
    Eigen::MatrixXd ErzRadial;

    // Integrating factor parameters 
    Eigen::MatrixXd GL,GR;
    Eigen::VectorXd g0,g1;

    // flux boundary conditions  
    Eigen::VectorXd wFluxBC;
    Eigen::VectorXd eFluxBC;
    Eigen::VectorXd nFluxBC;
    Eigen::VectorXd sFluxBC;

    // current boundary conditions
    Eigen::VectorXd wCurrentRBC;
    Eigen::VectorXd eCurrentRBC;
    Eigen::VectorXd nCurrentZBC;
    Eigen::VectorXd sCurrentZBC;

    // vectors for robin boundary conditions
    Eigen::VectorXd eInwardCurrentBC;
    Eigen::VectorXd nInwardCurrentBC;
    Eigen::VectorXd sInwardCurrentBC;

    Eigen::VectorXd eInwardFluxBC;
    Eigen::VectorXd nInwardFluxBC;
    Eigen::VectorXd sInwardFluxBC;

    Eigen::VectorXd eOutwardCurrToFluxRatioBC;
    Eigen::VectorXd nOutwardCurrToFluxRatioBC;
    Eigen::VectorXd sOutwardCurrToFluxRatioBC;
    
    Eigen::VectorXd eOutwardCurrToFluxRatioInwardWeightedBC;
    Eigen::VectorXd nOutwardCurrToFluxRatioInwardWeightedBC;
    Eigen::VectorXd sOutwardCurrToFluxRatioInwardWeightedBC;

    Eigen::VectorXd eAbsCurrentBC;
    Eigen::VectorXd nAbsCurrentBC;
    Eigen::VectorXd sAbsCurrentBC;

    bool useMGQDSources = false;

    MultiPhysicsCoupledQD * mpqd; 
    MultiGroupQD * mgqd;

    // FUNCTIONS

    shared_ptr<GreyGroupSolver> GGSolver;
    shared_ptr<GreyGroupSolverBase> GGSolverBase;

    GreyGroupQD(Materials * myMaterials,\
        Mesh * myMesh,\
        YAML::Node * myInput,\
        MultiPhysicsCoupledQD * myMPQD);

    void buildLinearSystem();
    void buildSteadyStateLinearSystem();
    void printBCParams();
    void printEddingtons();
    void assignMGQDPointer(MultiGroupQD * myMGQD);

    /* PETSc functions */

    // Steady state
    void buildSteadyStateLinearSystem_p();
    
    // Transient
    void buildLinearSystem_p();

  private:
    Materials * materials;
    Mesh * mesh; 
    YAML::Node * input;
};

//==============================================================================

#endif
