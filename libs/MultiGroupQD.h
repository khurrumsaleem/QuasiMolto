#ifndef MULTIGROUPQD_H
#define MULTIGROUPQD_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Materials.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std;
using namespace arma;

class SingleGroupQD; // forward declaration
class QDSolver; // forward declaration
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

  private:
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;

};

//==============================================================================

#endif
