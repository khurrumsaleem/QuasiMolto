// File: MultiGroupTransport.cpp
// Purpose: define a class that manipulates each single group transport object
// Date: October 28, 2019

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Materials.h"
#include "SingleGroupTransport.h"
#include "MultiGroupTransport.h"
#include "StartingAngle.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
//! MultiGroupTransport class object constructor

MultiGroupTransport::MultiGroupTransport(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput)
{
  materials = myMaterials;
  mesh = myMesh;
  input = myInput;
  cout << "number of energy groups: " <<materials->nGroups << endl;
  for (int iGroups = 0; iGroups < materials->nGroups; ++iGroups){
    shared_ptr<SingleGroupTransport> newSGT (new SingleGroupTransport(iGroups,\
      this,materials,mesh,input));
    SGTs.push_back(std::move(newSGT)); 
  }
  startAngleSolve = std::make_shared<StartingAngle>(mesh,materials,input);
  calcSources(); 
  solveStartAngles();
};

//==============================================================================

//==============================================================================
//! solveStartAngles loop over energy groups and call starting angle solver

void MultiGroupTransport::solveStartAngles()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->solveStartAngle();
  }
};

//==============================================================================

//==============================================================================
//! calcSources loop over SGTS and make call to calc sources in each

void MultiGroupTransport::calcSources()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->calcSource();
    cout << "Source for group "<< SGTs[iGroup]->energyGroup << endl;
    cout << SGTs[iGroup]->q << endl;
  }
};

//==============================================================================


