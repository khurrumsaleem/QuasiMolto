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

bool MultiGroupTransport::calcFluxes()
{
  
  typedef Eigen::Array<bool,Eigen::Dynamic,1> VectorXb;
  VectorXb converged(SGTs.size());
  Eigen::VectorXd residuals(SGTs.size());
  double epsilon = 1E-5;
  bool allConverged;

  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    residuals(iGroup)=SGTs[iGroup]->calcFlux();
    cout << "Scalar flux group " << iGroup << endl;
    cout << SGTs[iGroup]->sFlux<<endl;
    cout << endl;
    converged(iGroup) = residuals(iGroup) < epsilon;
  }

  allConverged = not(converged.isZero());

  return allConverged;  

};

//==============================================================================



//==============================================================================
//! calcSources loop over SGTS and make call to calc sources in each

void MultiGroupTransport::calcSources(string calcType)
{
  typedef Eigen::Array<bool,Eigen::Dynamic,1> VectorXb;
  VectorXb converged(SGTs.size());
  Eigen::VectorXd sourceNorms(SGTs.size());
  double epsilon = 1E-5;
  
  converged.fill(false); 
  sourceNorms.setZero();

  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->calcSource(calcType);
    //cout << "Source for group "<< SGTs[iGroup]->energyGroup << endl;
  }
};

//==============================================================================

//==============================================================================
//! calcSources loop over SGTS and make call to calc sources in each

void MultiGroupTransport::sourceIteration()
{
  typedef Eigen::Array<bool,Eigen::Dynamic,1> VectorXb;
  VectorXb converged(SGTs.size());
  cout << "In source iteration loop!"<< endl;
  calcSources(); 
  cout << "Fission source calculated!"<< endl;
  bool allConverged= false;
  int maxIter = 500;

  for (int iter = 0; iter < maxIter; ++iter){
  
    calcSources("s"); 
    solveStartAngles();
    solveSCBs();
    allConverged=calcFluxes();
    if (allConverged) {
      cout << "Converged in " << iter << " iterations."<< endl;
      break;
    }

  } 

  if(not(allConverged)){
    cout << "Fixed source iteration did NOT converge within " << maxIter;
    cout << " iterations." << endl;
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

