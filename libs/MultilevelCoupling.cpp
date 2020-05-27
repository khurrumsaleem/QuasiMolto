// File: MultilevelCoupling.cpp     
// Purpose: Solve a coupled MGT, MGQD, MPQD problem. 
// Date: May 27, 2020

#include "MultilevelCoupling.h"

using namespace std;

//==============================================================================
/// MultilevelCoupling class object constructor
///
/// @param [in] myMesh mesh object 
/// @param [in] myMats materials object 
/// @param [in] myInput input object 
/// @param [in] myMPQD MultiPhysicsCoupledQD object 
/// @param [in] myMGPQD multigroup quasidiffusion object 
MultilevelCoupling::MultilevelCoupling(Mesh * myMesh,\
        Materials * myMats,\
        YAML::Node * myInput,\
        MultiPhysicsCoupledQD * myMPQD,\
        MultiGroupQD * myMGQD)
{

  // Assign inputs to their member variables
  mesh = myMesh; 
  mats = myMats;
  input = myInput;
  mpqd = myMPQD;
  mgqd = myMGQD;

  // Assign MGQD pointer in grey group solver
  mgqd->assignMultiPhysicsCoupledQDPointer(mpqd);

  // Assign MGQD pointer in grey group solver
  mpqd->ggqd->assignMGQDPointer(mgqd);

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
bool MultilevelCoupling::solveOneStep()
{

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
void MultilevelCoupling::solveTransient()
{

};
//==============================================================================

