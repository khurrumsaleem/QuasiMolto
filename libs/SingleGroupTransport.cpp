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
  aHalfFlux.set_size(mesh->zCent.size(),mesh->rCent.size(),mesh->quadrature.size());
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
  cout << aFlux << endl;
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
      }
    }
  } 
};

//==============================================================================


