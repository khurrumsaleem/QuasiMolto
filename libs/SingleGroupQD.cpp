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
  // vectors for reading in optional parameters
  vector<double> inpaFlux, inpCurrent;
    
  // Assign energy group for this object 
  energyGroup = myEnergyGroup;
  
  // Assign pointers
  MGQD = myMGQD;
  mats = myMaterials;
  mesh = myMesh;
  input = myInput;

  // initialize Eddington factors to diffusion physics
  Err.setOnes(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  Err = (1.0/3.0)*Err;
  Ezz.setOnes(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  Ezz = (1.0/3.0)*Ezz;
  Erz.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  
  // initialize source 
  q.setOnes(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  q = 0.0*q;

  // initialize flux and current matrices
  sFlux.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  sFluxR.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size()+1);
  sFluxZ.setZero(mesh->zCornerCent.size()+1,mesh->rCornerCent.size());
  currentR.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size()+1);
  currentZ.setZero(mesh->zCornerCent.size()+1,mesh->rCornerCent.size());
  
  // initialize boundary conditions
  wFluxBC.setOnes(mesh->dzsCorner.size());
  eFluxBC.setOnes(mesh->dzsCorner.size());
  nFluxBC.setOnes(mesh->drsCorner.size());
  sFluxBC.setOnes(mesh->drsCorner.size());
  wCurrentRBC.setOnes(mesh->dzsCorner.size());
  eCurrentRBC.setOnes(mesh->dzsCorner.size());
  nCurrentZBC.setOnes(mesh->drsCorner.size());
  sCurrentZBC.setOnes(mesh->drsCorner.size());

  // check for boundary condition specified in input file

  if ((*input)["parameters"]["lowerBC"]){

    inpaFlux=(*input)["parameters"]["lowerBC"]\
      .as<vector<double> >();

    if (inpaFlux.size() == 1)
      nFluxBC = 4.0*inpaFlux[0]*nFluxBC;
    else
      nFluxBC = 4.0*inpaFlux[energyGroup]*nFluxBC;
  }

  if ((*input)["parameters"]["upperBC"]){

    inpaFlux=(*input)["parameters"]["upperBC"]\
      .as<vector<double> >();

    if (inpaFlux.size() == 1)
      sFluxBC = 4.0*inpaFlux[0]*sFluxBC;
    else
      sFluxBC = 4.0*inpaFlux[energyGroup]*sFluxBC;
  }

  if ((*input)["parameters"]["innerBC"]){

    inpaFlux=(*input)["parameters"]["innerBC"]\
      .as<vector<double> >();

    if (inpaFlux.size() == 1)
      wFluxBC = 4.0*inpaFlux[0]*wFluxBC;
    else
      wFluxBC = 4.0*inpaFlux[energyGroup]*wFluxBC;
  }

  if ((*input)["parameters"]["outerBC"]){

    inpaFlux=(*input)["parameters"]["outerBC"]\
      .as<vector<double> >();

    if (inpaFlux.size() == 1)
      eFluxBC = 4.0*inpaFlux[0]*eFluxBC;
    else
      eFluxBC = 4.0*inpaFlux[energyGroup]*eFluxBC;
  }

  if ((*input)["parameters"]["lowerCurrentBC"]){

    inpCurrent=(*input)["parameters"]["lowerCurrentBC"]\
      .as<vector<double> >();

    if (inpaFlux.size() == 1)
      nCurrentZBC = inpCurrent[0]*nCurrentZBC;
    else
      nCurrentZBC = inpCurrent[energyGroup]*nCurrentZBC;
  }

  if ((*input)["parameters"]["upperCurrentBC"]){

    inpCurrent=(*input)["parameters"]["upperCurrentBC"]\
      .as<vector<double> >();

    if (inpaFlux.size() == 1)
      sCurrentZBC = inpCurrent[0]*sCurrentZBC;
    else
      sCurrentZBC = inpCurrent[energyGroup]*sCurrentZBC;
  }

  if ((*input)["parameters"]["innerCurrentBC"]){

    inpCurrent=(*input)["parameters"]["innerCurrentBC"]\
      .as<vector<double> >();

    if (inpaFlux.size() == 1)
      wCurrentRBC = inpCurrent[0]*wCurrentRBC;
    else
      wCurrentRBC = inpCurrent[energyGroup]*wCurrentRBC;
  }

  if ((*input)["parameters"]["outerCurrentBC"]){

    inpCurrent=(*input)["parameters"]["outerCurrentBC"]\
      .as<vector<double> >();

    if (inpaFlux.size() == 1)
      eCurrentRBC = inpCurrent[0]*eCurrentRBC;
    else
      eCurrentRBC = inpCurrent[energyGroup]*eCurrentRBC;
  }



};
//==============================================================================

//==============================================================================
void SingleGroupQD::formContributionToLinearSystem()
{
  MGQD->QDSolve->formLinearSystem(this);
}
//==============================================================================

//==============================================================================
void SingleGroupQD::getFlux()
{
  MGQD->QDSolve->getFlux(this);
}
//==============================================================================

//==============================================================================
Eigen::VectorXd SingleGroupQD::getSolutionVector()
{
  return MGQD->QDSolve->getSolutionVector(this);
}
//==============================================================================
