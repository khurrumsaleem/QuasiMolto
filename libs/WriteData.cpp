/// File: WriteData.cpp
// Purpose: output data 
// Date: May 26, 2019

#include "WriteData.h"
#include "Mesh.h"

using namespace std; 

///==============================================================================
/// WriteData class object constructor
///
/// @param [in] nZ number of axial cells
/// @param [in] nR number of radial cells
WriteData::WriteData(Mesh * myMesh)
{
  // assign pointers
  mesh = myMesh;
};
//==============================================================================

