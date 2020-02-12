// File: MultiGroupQD.cpp
// Purpose: define a class that manipulates each single group quasidiffusion object
// Date: February 05, 2020

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Materials.h"
#include "SingleGroupQD.h"
#include "MultiGroupQD.h"
#include "QuasidiffusionSolver.h"
#include "MMS.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std;
using namespace arma;

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
void MultiGroupQD::buildLinearSystem()
{
  QDSolve->A.setZero();
  for (int iGroup = 0; iGroup < SGQDs.size(); iGroup++)
  {
    SGQDs[iGroup]->formContributionToLinearSystem();
  }
}
//==============================================================================

//==============================================================================
void MultiGroupQD::solveLinearSystem()
{
  QDSolve->solve();
  
  for (int iGroup = 0; iGroup < SGQDs.size(); iGroup++)
  {
    SGQDs[iGroup]->getFlux();
    cout << SGQDs[iGroup]->sFlux << endl;
  }
}
//==============================================================================

//==============================================================================
void MultiGroupQD::setInitialCondition()
{
  Eigen::VectorXd initialCondition(QDSolve->energyGroups*\
    QDSolve->nGroupUnknowns);
  initialCondition.setZero(); 
 
  for (int iGroup = 0; iGroup < SGQDs.size(); iGroup++)
  {
    initialCondition = initialCondition + SGQDs[iGroup]->getSolutionVector();
  }

  QDSolve->xPast = initialCondition;
}
//==============================================================================

//==============================================================================

void MultiGroupQD::solveMGQDOnly()
{
  setInitialCondition();
  for (int iTime = 0; iTime < mesh->dts.size(); iTime++)
  {
    buildLinearSystem();
    cout << "time: " <<mesh->ts[iTime+1] << endl;
    solveLinearSystem();
    QDSolve->xPast = QDSolve->x;
  }
  writeFluxes();
}
//==============================================================================

//==============================================================================
/// Wrapper over SGTs to write flux present in each

void MultiGroupQD::writeFluxes()
{
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup){
    SGQDs[iGroup]->writeFlux();
  }
};

//==============================================================================

