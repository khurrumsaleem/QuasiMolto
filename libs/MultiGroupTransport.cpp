// File: MultiGroupTransport.cpp
// Purpose: define a class that manipulates each single group transport object
// Date: October 28, 2019

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Material.h"
#include "SingleGroupTransport.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
//! MultiGroupTransport class that holds multigroup transport information

class MultiGroupTransport
{
  public:
  vector< shared_ptr<SingleGroupTransport> > SGTs;
  // public functions
  MultiGroupTransport(Mesh * myMesh,\
    YAML::Node * myInput);

  private:
  YAML::Node * input;
  Mesh * mesh;

};

//==============================================================================

//==============================================================================
//! MultiGroupTransport class object constructor

MultiGroupTransport::MultiGroupTransport(Mesh * myMesh,\
  YAML::Node * myInput)
{
  mesh = myMesh;
  input = myInput;
};

//==============================================================================

