#ifndef MULTIGROUPTRANSPORT_H
#define MULTIGROUPTRANSPORT_H

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

class SingleGroupTransport; // forward declaration
class StartingAngle; // forward declaration
class SimpleCornerBalance; // forward declaration

//==============================================================================
//! MultGroupTransport class that holds multigroup transport information

class MultiGroupTransport
{
  public:
  vector< shared_ptr<SingleGroupTransport> > SGTs;
  shared_ptr<StartingAngle> startAngleSolve;
  shared_ptr<SimpleCornerBalance> SCBSolve;
  // public functions
  MultiGroupTransport(Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput);
  void solveStartAngles();
  void calcSources();

  private:
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;
  
};

//==============================================================================


#endif
