// File: MultiGroupQDToMultiPhysicsQDCoupling.cpp     
// Purpose: Calculate coupling factors between multi-group quasidiffusion and
//   multi-physics quasidiffusion.
// Date: April 9, 2020

#include "MultiGroupQDToMultiPhysicsQDCoupling.h"

using namespace std;

//==============================================================================
/// MultiGroupQDToMultiPhysicsQDCoupling class object constructor
///
/// @param [in] blankType blank for this material
MGQDToMPQDCoupling::MGQDToMPQDCoupling(Mesh * myMesh,\
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

};
//==============================================================================
