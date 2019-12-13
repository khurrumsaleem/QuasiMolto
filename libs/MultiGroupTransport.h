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
  
  // default convergence criteria
  double epsAlpha=1E-3;
  double epsFlux=1E-5;
  double epsFissionSource=1E-5;
  
  // defaults for maximum iterations
  double sourceMaxIter=500;
  double powerMaxIter=10000;
  
  // public functions
  void solveStartAngles();
  void solveSCBs();
  bool calcSources(string calcType="FS");
  bool calcFluxes(string printResidual="noprint");
  bool calcAlphas(string printResidual="noprint");
  bool calcFissionSources(string printResidual="noprint");
  bool sourceIteration();
  bool powerIteration();
  void solveTransportOnly();
  void printDividers();
  void writeFluxes();

  private:
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;
  // for print formatting
  int spacing = 15;
  
};

//==============================================================================


#endif