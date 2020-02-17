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
  wFluxBC.setZero(mesh->dzsCorner.size());
  eFluxBC.setZero(mesh->dzsCorner.size());
  nFluxBC.setZero(mesh->drsCorner.size());
  sFluxBC.setZero(mesh->drsCorner.size());
  wCurrentRBC.setZero(mesh->dzsCorner.size());
  eCurrentRBC.setZero(mesh->dzsCorner.size());
  nCurrentZBC.setZero(mesh->drsCorner.size());
  sCurrentZBC.setZero(mesh->drsCorner.size());

  // check for optional parameters
  checkOptionalParams();
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
Eigen::VectorXd SingleGroupQD::getFluxSolutionVector()
{
  return MGQD->QDSolve->getFluxSolutionVector(this);
}
//==============================================================================

//==============================================================================
Eigen::VectorXd SingleGroupQD::getCurrentSolutionVector()
{
  return MGQD->QDSolve->getCurrentSolutionVector(this);
}
//==============================================================================


//==============================================================================
void SingleGroupQD::checkOptionalParams()
{
  // vectors for reading in optional parameters
  vector<double> inpaFlux,inpCurrent,inpSFluxPrev,inpCurrentPrev;

  // check for optional parameters specified in input file

  if ((*input)["parameters"]["lowerBC"]){

    inpaFlux=(*input)["parameters"]["lowerBC"]\
      .as<vector<double> >();
    
    nFluxBC.setOnes();
    if (inpaFlux.size() == 1)
      nFluxBC = 4.0*inpaFlux[0]*nFluxBC;
    else
      nFluxBC = 4.0*inpaFlux[energyGroup]*nFluxBC;
  }

  if ((*input)["parameters"]["upperBC"]){

    inpaFlux=(*input)["parameters"]["upperBC"]\
      .as<vector<double> >();

    sFluxBC.setOnes();
    if (inpaFlux.size() == 1)
      sFluxBC = 4.0*inpaFlux[0]*sFluxBC;
    else
      sFluxBC = 4.0*inpaFlux[energyGroup]*sFluxBC;
  }

  if ((*input)["parameters"]["innerBC"]){

    inpaFlux=(*input)["parameters"]["innerBC"]\
      .as<vector<double> >();

    wFluxBC.setOnes();
    if (inpaFlux.size() == 1)
      wFluxBC = 4.0*inpaFlux[0]*wFluxBC;
    else
      wFluxBC = 4.0*inpaFlux[energyGroup]*wFluxBC;
  }

  if ((*input)["parameters"]["outerBC"]){

    inpaFlux=(*input)["parameters"]["outerBC"]\
      .as<vector<double> >();

    eFluxBC.setOnes();
    if (inpaFlux.size() == 1)
      eFluxBC = 4.0*inpaFlux[0]*eFluxBC;
    else
      eFluxBC = 4.0*inpaFlux[energyGroup]*eFluxBC;
  }

  if ((*input)["parameters"]["lowerCurrentBC"]){

    inpCurrent=(*input)["parameters"]["lowerCurrentBC"]\
      .as<vector<double> >();

    nCurrentZBC.setOnes();
    if (inpaFlux.size() == 1)
      nCurrentZBC = inpCurrent[0]*nCurrentZBC;
    else
      nCurrentZBC = inpCurrent[energyGroup]*nCurrentZBC;
  }

  if ((*input)["parameters"]["upperCurrentBC"]){

    inpCurrent=(*input)["parameters"]["upperCurrentBC"]\
      .as<vector<double> >();

    sCurrentZBC.setOnes();
    if (inpaFlux.size() == 1)
      sCurrentZBC = inpCurrent[0]*sCurrentZBC;
    else
      sCurrentZBC = inpCurrent[energyGroup]*sCurrentZBC;
  }

  if ((*input)["parameters"]["innerCurrentBC"]){

    inpCurrent=(*input)["parameters"]["innerCurrentBC"]\
      .as<vector<double> >();

    wCurrentRBC.setOnes();
    if (inpaFlux.size() == 1)
      wCurrentRBC = inpCurrent[0]*wCurrentRBC;
    else
      wCurrentRBC = inpCurrent[energyGroup]*wCurrentRBC;
  }

  if ((*input)["parameters"]["outerCurrentBC"]){

    inpCurrent=(*input)["parameters"]["outerCurrentBC"]\
      .as<vector<double> >();

    eCurrentRBC.setOnes();
    if (inpaFlux.size() == 1)
      eCurrentRBC = inpCurrent[0]*eCurrentRBC;
    else
      eCurrentRBC = inpCurrent[energyGroup]*eCurrentRBC;
  }

  if ((*input)["parameters"]["initial previous flux"]){

    inpSFluxPrev=(*input)["parameters"]["initial previous flux"]\
      .as<vector<double> >();

    sFlux.setOnes();
    sFluxR.setOnes();
    sFluxZ.setOnes();
    if (inpSFluxPrev.size() == 1)
    {
      sFlux = inpSFluxPrev[0]*sFlux;
      sFluxR = inpSFluxPrev[0]*sFluxR;
      sFluxZ = inpSFluxPrev[0]*sFluxZ;
    }
    else
    {
      sFlux = inpSFluxPrev[energyGroup]*sFlux;
      sFluxR = inpSFluxPrev[energyGroup]*sFluxR;
      sFluxZ = inpSFluxPrev[energyGroup]*sFluxZ;
    }
  }

  if ((*input)["parameters"]["initial previous current"]){

    inpCurrentPrev=(*input)["parameters"]["initial previous current"]\
      .as<vector<double> >();

    currentR.setOnes();
    currentZ.setOnes();
    if (inpSFluxPrev.size() == 1)
    {
      currentR = inpCurrentPrev[0]*sFluxR;
      currentZ = inpCurrentPrev[0]*sFluxZ;
    }
    else
    {
      currentR = inpCurrentPrev[energyGroup]*sFluxR;
      currentZ = inpCurrentPrev[energyGroup]*sFluxZ;
    }
  }

}
//==============================================================================

//==============================================================================
/// Write the flux in this energy group to a CVS

void SingleGroupQD::writeFlux()
{

  ofstream fluxFile;
  string fileName;

  // parse file name
  fileName = "scalar-flux-group-" + to_string(energyGroup)\
    +"-" + to_string(mesh->dz)+ ".csv";

  // open file
  fluxFile.open(fileName);

  // write flux values to .csv
  for (int iZ = 0; iZ < sFlux.rows(); ++iZ) {
    fluxFile << sFlux(iZ,0);
    for (int iR = 1; iR < sFlux.cols(); ++iR) {
      fluxFile <<","<< sFlux(iZ,iR);
    }
    fluxFile << endl;
  }
  fluxFile.close();

  // if this is the first energy group, write mesh too
  if (energyGroup==0){
    // write radial mesh to .csv
    fileName = "r-mesh.csv";
    fluxFile.open(fileName);
    fluxFile << mesh->rCornerCent(0);
    for (int iR = 1; iR < mesh->rCornerCent.size(); ++iR) {
      fluxFile << "," << mesh->rCornerCent(iR);
    }
    fluxFile << endl;
    fluxFile.close();

    // write axial mesh to .csv
    fileName = "z-mesh.csv";
    fluxFile.open(fileName);
    fluxFile << mesh->zCornerCent(0);
    for (int iZ = 1; iZ < mesh->zCornerCent.size(); ++iZ) {
      fluxFile << "," << mesh->zCornerCent(iZ);
    }
    fluxFile << endl;
    fluxFile.close();
  }
};

//==============================================================================

