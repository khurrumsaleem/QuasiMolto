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
#include "SimpleCornerBalance.h"
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
  SCBSolve = std::make_shared<SimpleCornerBalance>(mesh,materials,input);
  sourceIteration();
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
//! solveSCBs loop over energy groups and call starting angle solver

void MultiGroupTransport::solveSCBs()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->solveSCB();
  }
};

//==============================================================================

//==============================================================================
//! solveSCBs loop over energy groups and call starting angle solver

void MultiGroupTransport::calcFluxes()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->calcFlux();
  }
};

//==============================================================================



//==============================================================================
//! calcSources loop over SGTS and make call to calc sources in each

void MultiGroupTransport::calcSources()
{
  typedef Eigen::Array<bool,Eigen::Dynamic,1> VectorXb;
  VectorXb converged(SGTs.size());
  Eigen::VectorXd sourceNorms(SGTs.size());
  double epsilon = 1E-5;

  converged.fill(false); 
  sourceNorms.setZero();

  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    sourceNorms(iGroup)=SGTs[iGroup]->calcSource();
    //cout << "Source for group "<< SGTs[iGroup]->energyGroup << endl;
    cout << SGTs[iGroup]->q << endl;
    converged(iGroup)=sourceNorms[iGroup]<epsilon;
  }
  cout << "sourceNorms: " << endl;
  cout << sourceNorms << endl;
};

//==============================================================================

//==============================================================================
//! calcSources loop over SGTS and make call to calc sources in each

void MultiGroupTransport::sourceIteration()
{
  typedef Eigen::Array<bool,Eigen::Dynamic,1> VectorXb;
  VectorXb converged(SGTs.size());

  for (int iter = 0; iter < 50; ++iter){
  
    calcSources(); 
    solveStartAngles();
    solveSCBs();
    calcFluxes();

  } 

  writeFluxes(); 
};

//==============================================================================

//==============================================================================
//! solveSCBs loop over energy groups and call starting angle solver

void MultiGroupTransport::writeFluxes()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->writeFlux();
  }
};

//==============================================================================

