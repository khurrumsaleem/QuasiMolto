// File: SingleGroupTransport.cpp
// Purpose: define a single group transport object
// Date: October 28, 2019

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Material.h"
#include "Materials.h"
#include "MultiGroupTransport.h"
#include "SingleGroupTransport.h"
#include "StartingAngle.h"
#include "SimpleCornerBalance.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
//! SingleGroupTransport class object constructor

SingleGroupTransport::SingleGroupTransport(int myEnergyGroup,\
  MultiGroupTransport * myMGT,\
  Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput)
{
  energyGroup = myEnergyGroup;
  MGT = myMGT;
  mats = myMaterials;
  mesh = myMesh;
  input = myInput;

  aFlux.set_size(mesh->zCent.size(),mesh->rCent.size(),mesh->nAngles);
  aFlux.zeros();
  aHalfFlux.set_size(mesh->zCent.size(),mesh->rCent.size(),\
    mesh->quadrature.size());
  aHalfFlux.zeros();
  sFlux.setOnes(mesh->zCent.size(),mesh->rCent.size());
  q.setZero(mesh->zCent.size(),mesh->rCent.size());
  cout << "Created transport energy group " << energyGroup << endl;
};

//==============================================================================

//==============================================================================
//! solveStartAngle call starting angle solver

void SingleGroupTransport::solveStartAngle()
{
  MGT->startAngleSolve->calcStartingAngle(&aHalfFlux,&q,energyGroup);
  cout << aHalfFlux << endl;
};

//==============================================================================

//==============================================================================
//! solveSCB call starting angle solver

void SingleGroupTransport::solveSCB()
{
  MGT->SCBSolve->solve(&aFlux,&aHalfFlux,&q,energyGroup);
  calcFlux();
  writeFlux();
};

//==============================================================================

//==============================================================================
//! calcSource calculate the source in an energy group

// Assuming the fluxes currently contained in each SGT object
void SingleGroupTransport::calcSource()
{
  q.setZero(mesh->zCent.size(),mesh->rCent.size());

  for (int iZ = 0; iZ < mesh->zCent.size(); ++iZ){

    for (int iR = 0; iR < mesh->rCent.size(); ++iR){

      for (int iGroup = 0; iGroup < MGT->SGTs.size(); ++iGroup){

        q(iZ,iR) = q(iZ,iR) \
        +mats->sigS(iZ,iR,iGroup,energyGroup)*MGT->SGTs[iGroup]->sFlux(iZ,iR)\
        + mats->nu(iZ,iR)*mats->sigF(iZ,iR,iGroup)\
        *MGT->SGTs[iGroup]->sFlux(iZ,iR);
        // need to account for precursors, too.
        // and allow input of xi
      } // iGroup
    } // iR
  } // iZ 
};

//==============================================================================

//==============================================================================
//! calcFlux calculate the flux in this energy group

// Assuming the fluxes currently contained in each SGT object
void SingleGroupTransport::calcFlux()
{
  
  double weight; 
  int weightIdx=3,angIdx;

  sFlux.setZero(mesh->zCent.size(),mesh->rCent.size());

  for (int iQ = 0; iQ < mesh->quadrature.size(); ++iQ){
    for (int iP = 0; iP < mesh->quadrature[iQ].nOrd; ++iP){

      weight = mesh->quadrature[iQ].quad[iP][weightIdx];
      angIdx = mesh->quadrature[iQ].ordIdx[iP];

      for (int iGroup = 0; iGroup < MGT->SGTs.size(); ++iGroup){
        for (int iZ = 0; iZ < mesh->zCent.size(); ++iZ){
          for (int iR = 0; iR < mesh->rCent.size(); ++iR){
            
            sFlux(iZ,iR) = sFlux(iZ,iR)\
            +weight*MGT->SGTs[iGroup]->aFlux(iZ,iR,angIdx);
          
          } // iR
        } // iZ
      } // iGroup
    } // iP
  } // iQ

  cout << "Scalar flux: " << endl;
  cout << sFlux << endl;
};

//==============================================================================

//==============================================================================
//! writeFlux write the flux in this energy group to a CVS

void SingleGroupTransport::writeFlux()
{

  ofstream fluxFile;
  string fileName;
  
  // parse file name 
  fileName = "scalar-flux-group-" + to_string(energyGroup) + ".csv";
  fluxFile.open(fileName);
  cout << "rows: " << sFlux.rows() << endl;
  cout << "cols: " << sFlux.cols() << endl;

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
    fluxFile << mesh->rCent(0);
    for (int iR = 1; iR < mesh->rCent.size(); ++iR) {
      fluxFile << "," << mesh->rCent(iR);
    }
    fluxFile << endl;
    fluxFile.close();
  
    // write axial mesh to .csv
    fileName = "z-mesh.csv"; 
    fluxFile.open(fileName); 
    fluxFile << mesh->zCent(0);
    for (int iZ = 1; iZ < mesh->zCent.size(); ++iZ) {
      fluxFile << "," << mesh->zCent(iZ);
    }
    fluxFile << endl;
    fluxFile.close();
  }
  
};

//==============================================================================




