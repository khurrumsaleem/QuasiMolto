// File: TransportToQDCoupling.cpp
// Purpose: handle communication between transport and quasidiffusion objects
// Date: February 12, 2020

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Materials.h"
#include "MultiGroupQD.h"
#include "MultiGroupTransport.h"
#include "TransportToQDCoupling.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std;
using namespace arma;

//==============================================================================
/// TransportToQDCoupling class object constructor
///
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myInput YAML input object for the simulation
TransportToQDCoupling::TransportToQDCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput,\
  MultiGroupTransport * myMGT,\
  MultiGroupQD * myMGQD)
{
  // Assign pointers for materials, mesh, and input objects
  materials = myMaterials;
  mesh = myMesh;
  input = myInput;

  // Assign pointers for multigroup transport and quasidiffusion objects
  MGT = myMGT;
  MGQD = myMGQD;

};
//==============================================================================
