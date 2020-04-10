// File: HeatTransfer.cpp     
// Purpose: Contains information and builds linear system for heat transfer
// Date: April 9, 2020

#include "HeatTransfer.h"

using namespace std;

//==============================================================================
/// HeatTransfer class object constructor
///
/// @param [in] blankType blank for this material
HeatTransfer::HeatTransfer(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput,\
  MultiPhysicsCoupledQD * myQD)
{
  // Assign inputs to their member variables
  mats = myMaterials;
  mesh = myMesh;
  input = myInput;
  qd = myQD;

  // Check for optional inputs 
  if ((*input)["parameters"]["wallTemp"]){
    wallT=(*input)["parameters"]["wallTemp"].as<double>();
  } 
  if ((*input)["parameters"]["inletTemp"]){
    inletT=(*input)["parameters"]["inletTemp"].as<double>();
  } 

  cout << "Initialized HeatTransfer object." << endl;
  cout << "wallTemp: " << wallT << endl;
  cout << "inletTemp: " << inletT << endl; 
};
//==============================================================================
