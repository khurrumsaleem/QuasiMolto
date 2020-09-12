// File: MultiGroupTransport.cpp
// Purpose: define a class that manipulates each single group transport object
// Date: October 28, 2019

#include "MultiGroupTransport.h"
#include "SingleGroupTransport.h"
#include "StartingAngle.h"
#include "SimpleCornerBalance.h"

using namespace std; 
using namespace arma;

//==============================================================================
/// MultiGroupTransport class object constructor
///
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myInput YAML input object for the simulation
MultiGroupTransport::MultiGroupTransport(Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput)
{
  // Assign pointers for materials, mesh, and input objects
  materials = myMaterials;
  mesh = myMesh;
  input = myInput;

  // Initialize pointers to each SGT group
  for (int iGroups = 0; iGroups < materials->nGroups; ++iGroups){
    shared_ptr<SingleGroupTransport> newSGT (new SingleGroupTransport(iGroups,\
          this,materials,mesh,input));
    SGTs.push_back(std::move(newSGT)); 
  }

  // Initialize StartingAngle and SimplCornerBalance solvers 
  startAngleSolve = std::make_shared<StartingAngle>(mesh,materials,input);
  SCBSolve = std::make_shared<SimpleCornerBalance>(mesh,materials,input);

  // Check to see if any convergence criteria are specified in input
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

};

//==============================================================================

//==============================================================================
/// Wrapper over SGTs to call starting angle solver

void MultiGroupTransport::solveStartAngles()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->solveStartAngle();
  }
};

//==============================================================================

//==============================================================================
/// Wrapper over SGTs to call SCB solver

void MultiGroupTransport::solveSCBs()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->solveSCB();
  }
};

//==============================================================================

//==============================================================================
/// Wrapper over SGTs to calculate scalar fluxes from angular fluxes
///
/// Uses the angular fluxes currently on each SGT object
/// @param [out] allConverged Indicates whether fluxes are globally converged
bool MultiGroupTransport::calcFluxes(string printResidual)
{

  // Vector containing the flux residuals in each SGT
  Eigen::VectorXd residuals(SGTs.size());
  residuals.setZero();

  // Vector indicating whether flux in each SGT is converged 
  VectorXb converged(SGTs.size());
  converged.fill(false); 

  // Criteria for flux convergence 
  double epsilon = 1E-5;

  // Boolean indicating whether all fluxes are converged
  bool allConverged=true;

  // Loop over SGTs, calculate fluxes, and determine whether the flux
  // in each SGT is converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    residuals(iGroup)=SGTs[iGroup]->calcFlux();
    converged(iGroup) = residuals(iGroup) < epsFlux;
  }

  // Print flux residuals
  if (printResidual == "print"){
    for (int iResidual = 0; iResidual < residuals.size(); ++iResidual){
      cout << setw(spacing) << residuals(iResidual);
    }
    cout << endl;
  }

  // Perform an AND operation over all SGTs to determine whether
  // flux is globally converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    allConverged = (allConverged and converged(iGroup)); 
  }

  return allConverged;  

};

//==============================================================================

//==============================================================================
/// Wrapper over SGTs to calculate sources 
///
/// @param [in] calcType Determines how source is calculated. "s" only 
/// re-evaluate scattering term. "fs" (default) re-evaluate scattering 
/// and fission terms
/// @param [out] allConverged Indicates whether sources are globally converged
bool MultiGroupTransport::calcSources(string calcType)
{

  // Vector indicating whether source in each SGT is converged
  VectorXb converged(SGTs.size());
  converged.fill(false); 

  // Vector containing the source residuals in each SGT
  Eigen::VectorXd residuals(SGTs.size());
  residuals.setZero();

  // Boolean indicating whether all sources are converged
  bool allConverged=true;

  // Loop over SGTs, calculate sources, and determine whether the source
  // in each SGT is converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    residuals(iGroup)=SGTs[iGroup]->calcSource(calcType);
    converged(iGroup) = residuals(iGroup) < epsFissionSource;
  }

  // Perform an AND operation over all SGTs to determine whether
  // sources are globally converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    allConverged = (allConverged and converged(iGroup)); 
  }

  return allConverged;
};

//==============================================================================

//==============================================================================
/// Wrapper over SGTs to calculate alphas
///
/// @param [out] allConverged Indicates whether alpha estimates are globally 
/// converged
bool MultiGroupTransport::calcAlphas(string printResidual,string calcType)
{

  // Vector indicating whether alpha in each SGT is converged
  VectorXb converged(SGTs.size());

  // Vector containing the alpha residuals in each SGT
  Eigen::VectorXd residuals(SGTs.size());

  // Boolean indicating whether all alphas are converged
  bool allConverged=true;

  // Loop over SGTs, calculate alphas, and determine whether the alphas
  // in each SGT are converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    residuals(iGroup)=SGTs[iGroup]->calcAlpha(calcType);
    converged(iGroup) = residuals(iGroup) < epsAlpha;
  }

  // Print alpha residuals
  if (printResidual == "print"){
    cout << "Alpha residuals: " << endl;
    for (int iResidual = 0; iResidual < residuals.size(); ++iResidual){
      cout << setw(spacing) << residuals(iResidual);
    }
    cout << endl; 
    cout << endl;
  } 

  // Perform an AND operation over all SGTs to determine whether
  // the alpha estimates are globally converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    allConverged = (allConverged and converged(iGroup)); 
  }

  return allConverged;  

};

//==============================================================================

//==============================================================================
/// Wrapper over SGTs to calculate the fission source.
///
/// @param [out] allConverged Indicates whether the newly calculated fission 
/// source is converged based on the convergence criteria epsilon  
bool MultiGroupTransport::calcFissionSources(string printResidual)
{

  // Vector indicating whether fission source in each SGT is converged
  VectorXb converged(SGTs.size());

  // Vector containing the fission source residuals in each SGT
  Eigen::VectorXd residuals(SGTs.size());

  // Criteria for fission source convergence 
  double epsilon = 1E-5;

  // Boolean indicating whether all fission sources are converged
  bool allConverged=true;

  // Loop over SGTs, calculate fission source, and determine whether the fission
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

  // Perform an AND operation over all SGTs to determine whether
  // fission sources are globally converged
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    allConverged = (allConverged and converged(iGroup)); 
  }

  return allConverged;  

};

//==============================================================================

//==============================================================================
/// Iterate on a solution using a fixed source
///
/// @param [out] allConverged Indicates whether fixed source iteration is
/// converged
bool MultiGroupTransport::sourceIteration()
{
  // Boolean indicated whether flux is globally converged
  bool allConverged= false;

  // Iterate on fission source
  for (int iter = 0; iter < sourceMaxIter; ++iter){

    // Calculate scatter source, solve for the starting angle equation, then 
    // solve the full time dependent neutron tranport equation with a 
    // fixed source. Then calculate the scalar flux with the newly calculated
    // angular flux
    calcSources("s"); 
    solveStartAngles();
    solveSCBs();
    allConverged=calcFluxes("print");

    // If the fluxes are globally converged, break out of the for loop
    if (allConverged) {
      cout << "Converged in " << iter+1 << " iteration(s)."<< endl;
      break;
    }

  } 

  // Print a statement indicating fixed source iteration was unsuccessful
  if(not(allConverged)){
    cout << "Fixed source iteration did NOT converge within " << sourceMaxIter;
    cout << " iterations." << endl;
  }

  return allConverged;
};

//==============================================================================

//==============================================================================
/// Perform a power iteration 
///
/// Update alphas and fission source in each SGT
/// @param [out] allConverged Indicates whether alpha and fission source 
/// estimates are converged
bool MultiGroupTransport::powerIteration()
{
  // Booleans indicating whether the alphas and fission source estimates are 
  // globally converged
  bool alphaConverged=false,fissionConverged=false,allConverged=false; 


  // Calculate alphas in each SGT
  alphaConverged = calcAlphas("print"); 

  // Calculate fission source in SGT
  fissionConverged=calcFissionSources("print");

  allConverged = (alphaConverged and fissionConverged);

  return allConverged;
};

//==============================================================================

//==============================================================================
/// Solve without any multiphysics coupling
///
/// Iterates using power (outer) and source (inner) iteration functions
void MultiGroupTransport::solveTransportOnly()
{
  // Booleans indicating whether source and power iterations are converged
  bool fluxConverged,fissionSourceConverged;

  // Maximum number of iterations
  int maxIter = 10000;

  // Calculate initial source
  calcSources(); 

  for (int iTime = 0; iTime < mesh->dts.size(); iTime++)
  {
    cout << "Flux residuals: " << endl;
    for (int iter = 0; iter < powerMaxIter; ++iter){

      // Perform source iteration
      fluxConverged=sourceIteration();

      if (fluxConverged){

        cout << endl; 

        // Perform power iteration
        printDividers();
        fissionSourceConverged=powerIteration(); 
        printDividers(); 

        // Problem is fully converged 
        if (fissionSourceConverged)
        {
          cout << "Solution converged!" << endl;
          break;
        } else
        {
          cout << "Flux residuals: " << endl;
        }

      } else 
      {
        cout << "Source iteration non-convergent." << endl;
        break;
      }
    }

    // set previous fluxes in each group
    for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup)
    {
      SGTs[iGroup]->sFluxPrev = SGTs[iGroup]->sFlux;
    }
  }
  writeFluxes(); 
};

//==============================================================================

//==============================================================================
/// Wrapper over SGTs to write flux present in each

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
/// Wrapper over SGTs to write flux present in each

void MultiGroupTransport::writeFluxes()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGTs[iGroup]->writeFlux();
  }
};

//==============================================================================

//==============================================================================
/// Extracts fluxes and currents from solution vector into 2D matrices 
void MultiGroupTransport::assignMultiPhysicsCoupledQDPointers\
       (MultiGroupQD * myMGQD, MultiPhysicsCoupledQD * myMPQD)
{
  mgqd = myMGQD;
  mpqd = myMPQD;
  useMPQDSources = true;
}
//==============================================================================


