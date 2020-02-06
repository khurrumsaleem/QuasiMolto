// File: QuasidiffusionSolver.cpp
// Purpose: Solve RZ quasidiffusion equations  
// Date: February 05, 2020

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Materials.h"
#include "Material.h"
#include "MMS.h"
#include "QuasidiffusionSolver.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
/// QuasidiffusionSolver object constructor
///
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myInput YAML input object for the simulation
QDSolver::QDSolver(Mesh * myMesh,\
  Materials * myMaterials,\
  YAML::Node * myInput)	      
{
  // Point to variables for mesh and input file
  mesh = myMesh;
  input = myInput;
  materials = myMaterials;

};

//==============================================================================

