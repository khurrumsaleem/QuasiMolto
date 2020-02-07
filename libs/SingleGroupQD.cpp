// File: SingleGroupQD.cpp
// Purpose: define a single group quasidiffusion object
// Date: February 5, 2020

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Material.h"
#include "Materials.h"
#include "MultiGroupQD.h"
#include "SingleGroupQD.h"
#include "QuasidiffusionSolver.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
/// Class object constructor
///
/// @param [in] myEnergyGroup Energy group for this single group transport
/// object
/// @param [in] myMGT MultiGroupTransport object for this simulation
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myInput YAML input object for the simulation
SingleGroupQD::SingleGroupQD(int myEnergyGroup,\
  MultiGroupQD * myMGQD,\
  Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput)
{
  // Assign energy group for this object 
  energyGroup = myEnergyGroup;
  
  // Assign pointers
  MGQD = myMGQD;
  mats = myMaterials;
  mesh = myMesh;
  input = myInput;

  // initialize Eddington factors to diffusion physics
  Err.setOnes(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  Err = (1/3)*Err;
  Ezz.setOnes(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  Ezz = (1/3)*Ezz;
  Erz.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size());

};
//==============================================================================

//==============================================================================
void SingleGroupQD::formContributionToLinearSystem()
{
  MGQD->QDSolve->formLinearSystem(this);
}
//==============================================================================
