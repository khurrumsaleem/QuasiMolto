// File: MultiGroupTransport.cpp
// Purpose: define a class that manipulates each single group transport object
// Date: October 28, 2019

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Materials.h"
#include "SingleGroupTransport.h"
#include "MultiGroupTransport.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
//! MultiGroupTransport class object constructor

MultiGroupTransport::MultiGroupTransport(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput)
{
  materials = myMaterials;
  mesh = myMesh;
  input = myInput;
 // cout << "number of energy groups: " <<materials->nGroups << endl;
};

//==============================================================================

