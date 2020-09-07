#ifndef MULTIGROUPQD_H
#define MULTIGROUPQD_H

#include "Mesh.h"
#include "QuasidiffusionSolver.h"
#include "MultiPhysicsCoupledQD.h"

// ToDo: this class should only refer to POINTERS of the forward declared type,
// but that is not the case. To work around in, the SingleGroupQD.h is included
// in the MultiGroupQD.cpp. Need to clean this up. Similar kluges are present
// in QuasidiffusionSolver, MultiGroupTransport, and SimpleCornerBalance.
class SingleGroupQD; // forward declaration

using namespace std;

//==============================================================================
//! MultGroupTransport class that holds multigroup transport information

class MultiGroupQD
{
  public:
    shared_ptr<QDSolver> QDSolve;
    vector< shared_ptr<SingleGroupQD> > SGQDs;
    // public functions
    MultiGroupQD(Materials * myMaterials,\
        Mesh * myMesh,\
        YAML::Node * myInput);
    void buildLinearSystem();
    void buildSteadyStateLinearSystem();
    void solveLinearSystem();
    void solveLinearSystemIterative();
    void buildBackCalcSystem();
    void buildSteadyStateBackCalcSystem();
    void backCalculateCurrent();
    void setInitialCondition();
    void solveMGQDOnly();
    void getFluxes();
    void updateVarsAfterConvergence();
    void updateSteadyStateVarsAfterConvergence();
    void writeFluxes();
    void writeVars();
    void printVars();
    void printEddingtons();
    void assignMultiPhysicsCoupledQDPointer(MultiPhysicsCoupledQD * myMPQD);
    string outputDir = "MGQD/";

  private:
    YAML::Node * input;
    Mesh * mesh;
    Materials * materials;

};

//==============================================================================

#endif
