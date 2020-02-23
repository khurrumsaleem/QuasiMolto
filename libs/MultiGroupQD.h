#ifndef MULTIGROUPQD_H
#define MULTIGROUPQD_H

#include <iomanip>
#include "Mesh.h"
#include "SingleGroupQD.h"
#include "QuasidiffusionSolver.h"

using namespace std;

//class SingleGroupQD; // forward declaration
//class QDSolver; // forward declaration
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
  void solveLinearSystem();
  void buildBackCalcSystem();
  void backCalculateCurrent();
  void setInitialCondition();
  void solveMGQDOnly();
  void writeFluxes();

  private:
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;

};

//==============================================================================

#endif
