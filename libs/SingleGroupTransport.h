#ifndef SINGLEGROUPTRANSPORT_H
#define SINGLEGROUPTRANSPORT_H

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

class MultiGroupTransport; // forward declaration

//==============================================================================
//! SingleGroupTransport class that holds multigroup transport information

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


#endif
