#ifndef MULTIGROUPTRANSPORT_H
#define MULTIGROUPTRANSPORT_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Material.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

class SingleGroupTransport; // forward declaration

//==============================================================================
//! MultGroupTransport class that holds multigroup transport information

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


#endif
