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
  cout << "Created transport energy group " << energyGroup << endl;
};

//==============================================================================

//==============================================================================
//! SingleGroupTransport class object constructor

void SingleGroupTransport::solveStartAngle()
{
  MGT->startAngleSolve->calcStartingAngle(&aHalfFlux,&sFlux,energyGroup);
};

//==============================================================================

