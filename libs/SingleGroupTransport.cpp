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
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
//! SingleGroupTransport class that holds transport information for one group

class SingleGroupTransport
{
  public:
  // public functions
  SingleGroupTransport(MultiGroupTransport * myMGT,\
    Mesh * myMesh,\
    YAML::Node * myInput);

  private:
  MultiGroupTransport * MGT;
  YAML::Node * input;
  Mesh * mesh;

};

//==============================================================================

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

