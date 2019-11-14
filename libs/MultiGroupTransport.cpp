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
  // assign pointers for materials, mesh, and input objects
  materials = myMaterials;
  mesh = myMesh;
  input = myInput;

  // initialize pointers to each SGT group
  for (int iGroups = 0; iGroups < materials->nGroups; ++iGroups){
    shared_ptr<SingleGroupTransport> newSGT (new SingleGroupTransport(iGroups,\
      this,materials,mesh,input));
    SGTs.push_back(std::move(newSGT)); 
  }

  // initialize StartingAngle and SimplCornerBalance solvers 
  startAngleSolve = std::make_shared<StartingAngle>(mesh,materials,input);
  SCBSolve = std::make_shared<SimpleCornerBalance>(mesh,materials,input);

  // check to see if any convergence criteria are specified in input
  if ((*input)["parameters"]["epsAlpha"]){
    epsAlpha=(*input)["parameters"]["epsAlpha"].as<double>();
  }
  if ((*input)["parameters"]["epsFlux"]){
    epsFlux=(*input)["parameters"]["epsFlux"].as<double>();
  }
  if ((*input)["parameters"]["epsFissionSource"]){
    epsFissionSource=(*input)["parameters"]["epsFissionSource"].as<double>();
  }
  if ((*input)["parameters"]["sourceMaxIter"]){
    epsAlpha=(*input)["parameters"]["sourceMaxIter"].as<double>();
  }
  if ((*input)["parameters"]["powerMaxIter"]){
    epsAlpha=(*input)["parameters"]["powerMaxIter"].as<double>();
  }
  
  // call source/power iteration solver
  solveTransportOnly();
};

//==============================================================================

//==============================================================================
//! solveStartAngles wrapper over SGTs to call starting angle solver

void MultiGroupTransport::solveStartAngles()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->solveStartAngle();
  }
};

//==============================================================================

//==============================================================================
//! solveSCBs wrapper over SGTs to call SCB solver

void MultiGroupTransport::solveSCBs()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->solveSCB();
  }
};

//==============================================================================

//==============================================================================
//! calcFluxes wrapper over SGTs to calculate scalar fluxes from angular fluxes

// uses the angular fluxes currently on each SGT object
bool MultiGroupTransport::calcFluxes(string printResidual)
{
  
  //typedef Eigen::Array<bool,Eigen::Dynamic,1> VectorXb;
  
  // vector containing the flux residuals in each SGT
  Eigen::VectorXd residuals(SGTs.size());
  residuals.setZero();
  
  // vector indicating whether flux in each SGT is converged 
  VectorXb converged(SGTs.size());
  converged.fill(false); 

  // criteria for flux convergence 
  double epsilon = 1E-5;

  // boolean indicating whether all fluxes are converged
  bool allConverged=true;

  // loop over SGTs, calculate fluxes, and determine whether the flux
  // in each SGT is converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    residuals(iGroup)=SGTs[iGroup]->calcFlux();
    converged(iGroup) = residuals(iGroup) < epsFlux;
  }

  // print flux residuals
  if (printResidual == "print"){
    for (int iResidual = 0; iResidual < residuals.size(); ++iResidual){
      cout << setw(spacing) << residuals(iResidual);
    }
    cout << endl;
  }
  //cout << setw(10) <<residuals.transpose() << endl;

  // perform an AND operation over all SGTs to determine whether
  // flux is globally converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    allConverged = (allConverged and converged(iGroup)); 
  }

  return allConverged;  

};

//==============================================================================

//==============================================================================
//! calcSources wrapper over SGTs to calculate sources 

// calcType determines how source is calculated
//   "s" only re-evaluate scattering term
//   "fs" (default) re-evaluate scattering and fission terms
bool MultiGroupTransport::calcSources(string calcType)
{
  
  // vector indicating whether source in each SGT is converged
  VectorXb converged(SGTs.size());
  converged.fill(false); 
  
  // vector containing the source residuals in each SGT
  Eigen::VectorXd residuals(SGTs.size());
  residuals.setZero();

  // criteria for source convergence 
  double epsilon = 1E-5;
  
  // boolean indicating whether all sources are converged
  bool allConverged=true;

  // loop over SGTs, calculate sources, and determine whether the source
  // in each SGT is converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    residuals(iGroup)=SGTs[iGroup]->calcSource(calcType);
    converged(iGroup) = residuals(iGroup) < epsFissionSource;
  }

  // perform an AND operation over all SGTs to determine whether
  // sources are globally converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    allConverged = (allConverged and converged(iGroup)); 
  }
  
  return allConverged;
};

//==============================================================================

//==============================================================================
//! calcAlphas wrapper over SGTs to calculate alphas

// In line with the alpha approximation to the time dependent neutron transport 
// equation
bool MultiGroupTransport::calcAlphas(string printResidual)
{
  
  // vector indicating whether alpha in each SGT is converged
  VectorXb converged(SGTs.size());

  // vector containing the alpha residuals in each SGT
  Eigen::VectorXd residuals(SGTs.size());

  // criteria for alpha convergence 
  double epsilon = 1E-3;

  // boolean indicating whether all alphas are converged
  bool allConverged=true;

  // loop over SGTs, calculate alphas, and determine whether the alphas
  // in each SGT are converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    residuals(iGroup)=SGTs[iGroup]->calcAlpha();
    converged(iGroup) = residuals(iGroup) < epsAlpha;
  }

  // print alpha residuals
  if (printResidual == "print"){
    cout << "Alpha residuals: " << endl;
    for (int iResidual = 0; iResidual < residuals.size(); ++iResidual){
      cout << setw(spacing) << residuals(iResidual);
    }
    cout << endl; 
    cout << endl;
  } 
 
  // perform an AND operation over all SGTs to determine whether
  // the alpha estimates are globally converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    allConverged = (allConverged and converged(iGroup)); 
  }

  return allConverged;  

};

//==============================================================================

//==============================================================================
//! calcFissionSources wrapper over SGTs to calculate the fission source.

// allConverged: indicates whether the newly calculated fission source is 
//   converged based on the convergence criteria epsilon  
bool MultiGroupTransport::calcFissionSources(string printResidual)
{
  
  // vector indicating whether fission source in each SGT is converged
  VectorXb converged(SGTs.size());

  // vector containing the fission source residuals in each SGT
  Eigen::VectorXd residuals(SGTs.size());

  // criteria for fission source convergence 
  double epsilon = 1E-5;

  // boolean indicating whether all fission sources are converged
  bool allConverged=true;

  // loop over SGTs, calculate fission source, and determine whether the fission
  // sources in each SGT are converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    residuals(iGroup)=SGTs[iGroup]->calcFissionSource();
    converged(iGroup) = residuals(iGroup) < epsFissionSource;
  }

  if (printResidual == "print"){
    cout << "Fission source residuals: " << endl;
    for (int iResidual = 0; iResidual < residuals.size(); ++iResidual){
      cout << setw(spacing) << residuals(iResidual);
    }
    cout << endl; 
    cout << endl; 
  }

  // perform an AND operation over all SGTs to determine whether
  // fission sources are globally converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    allConverged = (allConverged and converged(iGroup)); 
  }

  return allConverged;  

};

//==============================================================================

//==============================================================================
//! sourceIteration iterate on a solution for a fixed source

// allConverged: indicates whether fluxes are converged
bool MultiGroupTransport::sourceIteration()
{
  // boolean indicated whether flux is globally converged
  bool allConverged= false;

  // iterate on fission source
  for (int iter = 0; iter < sourceMaxIter; ++iter){
    
    // calculate scatter source, solve for the starting angle equation, then 
    // solve the full time dependent neutron tranport equation with a 
    // fixed source. Then calculate the scalar flux with the newly calculated
    // angular flux
    calcSources("s"); 
    solveStartAngles();
    solveSCBs();
    allConverged=calcFluxes("print");

    // if the fluxes are globally converged, break out of the for loop
    if (allConverged) {
      cout << "Converged in " << iter << " iterations."<< endl;
      break;
    }

  } 

  // print a statement indicating fixed source iteration was unsuccessful
  if(not(allConverged)){
    cout << "Fixed source iteration did NOT converge within " << sourceMaxIter;
    cout << " iterations." << endl;
  }

  return allConverged;
};

//==============================================================================

//==============================================================================
//! powerIteration perform a power iteration 

// converged: indicates whether alpha and fission source estimates are converged
// update alphas and fission source in each SGT
bool MultiGroupTransport::powerIteration()
{
  // booleans indicating whether the alphas and fission source estimates are 
  // globally converged
  bool alphaConverged=false,fissionConverged=false; 

  // calculate alphas in each SGT
  alphaConverged=calcAlphas("print");

  // calculate fission source in SGT
  fissionConverged=calcFissionSources("print"); 

  return (alphaConverged and fissionConverged);
};

//==============================================================================

//==============================================================================
//! solveTransportOnly solve without any multiphysics coupling

// iterates using power (outer) and source (inner) iteration functions
void MultiGroupTransport::solveTransportOnly()
{
  // booleans indicating whether source and power iterations are converged
  bool fluxConverged,fissionSourceConverged;

  // maximum number of iterations
  int maxIter = 10000;

  // calculate initial source
  calcSources(); 

  cout << "Flux residuals: " << endl;
  for (int iter = 0; iter < powerMaxIter; ++iter){

    // perform source iteration
    fluxConverged=sourceIteration();

    if (fluxConverged){

      cout << endl; 

      // perform power iteration
      printDividers();
      fissionSourceConverged=powerIteration(); 
      printDividers(); 

      // problem is fully converged 
      if (fissionSourceConverged){
        cout << "Solution converged!" << endl;
        break;
      } else{
        cout << "Flux residuals: " << endl;
      }
    
    } else {
      cout << "Source iteration non-convergent." << endl;
      break;
    }
  } 

  writeFluxes(); 
};

//==============================================================================

//==============================================================================
//! writeFluxes wrapper over SGTs to write flux present in each

void MultiGroupTransport::printDividers()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    for (int iSpace = 0; iSpace < spacing; ++iSpace ){
      cout << "=";
    }
  }
  cout << endl;
  cout << endl;
};

//==============================================================================


//==============================================================================
//! writeFluxes wrapper over SGTs to write flux present in each

void MultiGroupTransport::writeFluxes()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->writeFlux();
  }
};

//==============================================================================

