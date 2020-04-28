// File: GreyGroupQD.cpp     
// Purpose: Contains information and builds linear system for grey group qd
// Date: April 9, 2020

#include "GreyGroupQD.h"

using namespace std;

//==============================================================================
/// GreyGroupQD class object constructor
///
/// @param [in] blankType blank for this material
GreyGroupQD::GreyGroupQD(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput)
{
  // Assign pointers for materials, mesh, and input objects
  materials = myMaterials;
  mesh = myMesh;
  input = myInput;

  GGSolver = std::make_shared<GreyGroupSolver>(this,mesh,materials,input);

};
//==============================================================================
