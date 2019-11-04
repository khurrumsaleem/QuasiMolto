// File: SingleGroupTransport.cpp
// Purpose: define a single group transport object
// Date: October 28, 2019

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Material.h"
#include "MultiGroupTransport.h"
#include "SingleGroupTransport.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
//! SingleGroupTransport class object constructor

SingleGroupTransport::SingleGroupTransport(MultiGroupTransport * myMGT,\
  Mesh * myMesh,\
  YAML::Node * myInput)
{
  MGT = myMGT;
  mesh = myMesh;
  input = myInput;
};

//==============================================================================

