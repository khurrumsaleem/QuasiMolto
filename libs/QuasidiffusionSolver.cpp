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

  // temporary variables for initialization
  int nUnknowns;

  // calculate number of unknowns  
  energyGroups = materials->nGroups;
  nUnknowns = energyGroups*(5*(mesh->zCornerCent.size()\
    *mesh->rCornerCent.size()) + 2*mesh->zCornerCent.size()\
    + 2*mesh->rCornerCent.size());

  // initialize size of linear system
  A.resize(nUnknowns,nUnknowns);

};

//==============================================================================

//==============================================================================
/// Form a portion of the linear system that belongs to SGQD 
///
void QDSolver::formLinearSystem(SingleGroupQD * SGQD)	      
{

  cout << "linear system formed" << endl;

};

//==============================================================================


