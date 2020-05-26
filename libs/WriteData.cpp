/// File: WriteData.cpp
// Purpose: output data 
// Date: May 26, 2019

#include "Mesh.h"
#include "WriteData.h"

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
//=============================================================================

///==============================================================================
/// WriteData class object constructor
///
/// @param [in] dirName relative path of directory to be made 
void WriteData::makeDirectory(string myDirName)
{
 
  // Get time at present state
  string timeDir = outputDirectory + to_string(mesh->ts[mesh->state]) + "/";
 
  // Parse system command
  string command = "mkdir -p " + timeDir + myDirName;  
 
  // Issue system command to create directory
  system(command.c_str());

};
//==============================================================================

