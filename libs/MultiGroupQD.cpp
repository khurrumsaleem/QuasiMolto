// File: MultiGroupQD.cpp
// Purpose: define a class that manipulates each single group quasidiffusion object
// Date: February 05, 2020

#include "MultiGroupQD.h"
#include "SingleGroupQD.h"

using namespace std;

//==============================================================================
/// MultiGroupQD class object constructor
///
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myInput YAML input object for the simulation
MultiGroupQD::MultiGroupQD(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput)
{
  // Assign pointers for materials, mesh, and input objects
  materials = myMaterials;
  mesh = myMesh;
  input = myInput;

  // Initialize pointers to each SGQD group
  for (int iGroups = 0; iGroups < materials->nGroups; ++iGroups){
    shared_ptr<SingleGroupQD> newSGQD (new SingleGroupQD(iGroups,\
      this,materials,mesh,input));
    SGQDs.push_back(std::move(newSGQD));
  }

  QDSolve = std::make_shared<QDSolver>(mesh,materials,input);

};
//==============================================================================

//==============================================================================
/// Loops over energy groups and builds the linear system to solve the 
/// multigroup quasidiffusion equations
void MultiGroupQD::buildLinearSystem()
{
  QDSolve->A.setZero();
  QDSolve->b.setZero();
  for (int iGroup = 0; iGroup < SGQDs.size(); iGroup++)
  {
    SGQDs[iGroup]->formContributionToLinearSystem();
  }
}
//==============================================================================

//==============================================================================
/// Solves the linear system formed by the muligroup quasidiffusion equations
void MultiGroupQD::solveLinearSystem()
{
  QDSolve->solve();
  
  for (int iGroup = 0; iGroup < SGQDs.size(); iGroup++)
  {
    SGQDs[iGroup]->getFlux();
   // cout << SGQDs[iGroup]->sFlux << endl;
  }
}
//==============================================================================

//==============================================================================
/// Loops over energy groups and builds linear system to calculate the net 
/// neutron currents from the flux values currrently held in x, the solution
/// vector
void MultiGroupQD::buildBackCalcSystem()
{
  QDSolve->C.setZero();
  QDSolve->d.setZero();
  for (int iGroup = 0; iGroup < SGQDs.size(); iGroup++)
  {
    SGQDs[iGroup]->formContributionToBackCalcSystem();
  }
}
//==============================================================================

//==============================================================================
/// Solve the linear system which calculates net currents from the current flux
/// values
void MultiGroupQD::backCalculateCurrent()
{
  QDSolve->backCalculateCurrent();
  
//  for (int iGroup = 0; iGroup < SGQDs.size(); iGroup++)
//  {
//    SGQDs[iGroup]->getFlux();
//    cout << SGQDs[iGroup]->sFlux << endl;
//  }
}
//==============================================================================


//==============================================================================
/// Set the initial previous solution vectors to the values currently held in 
/// the flux and current matrices
void MultiGroupQD::setInitialCondition()
{
  // Initialize vectors
  Eigen::VectorXd initialFluxCondition(QDSolve->energyGroups*\
    QDSolve->nGroupUnknowns);
  Eigen::VectorXd initialCurrentCondition(QDSolve->energyGroups*\
    QDSolve->nGroupCurrentUnknowns);
  initialFluxCondition.setZero();   
  initialCurrentCondition.setZero(); 

  // Populate initial vectors with values from each energy group  
  for (int iGroup = 0; iGroup < SGQDs.size(); iGroup++)
  {
    initialFluxCondition = initialFluxCondition\
      + SGQDs[iGroup]->getFluxSolutionVector();
    initialCurrentCondition = initialCurrentCondition\
      + SGQDs[iGroup]->getCurrentSolutionVector();
  }

  // Set initial conditions in QDSolver object
  QDSolve->xPast = initialFluxCondition;
  QDSolve->currPast = initialCurrentCondition;
}
//==============================================================================

//==============================================================================
/// Solve a transient problem without any transport coupling using diffusion
/// values for the Eddington factors
void MultiGroupQD::solveMGQDOnly()
{
  setInitialCondition();
  for (int iTime = 0; iTime < mesh->dts.size(); iTime++)
  {
    buildLinearSystem();
    cout << "time: " <<mesh->ts[iTime+1] << endl;
    solveLinearSystem();
    QDSolve->xPast = QDSolve->x;
    buildBackCalcSystem();
    backCalculateCurrent();
  }
  writeFluxes();
}
//==============================================================================

//==============================================================================
/// Wrapper over SGQDs to write flux present in each
void MultiGroupQD::writeFluxes()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGQDs[iGroup]->writeFlux();
  }
}
//==============================================================================

