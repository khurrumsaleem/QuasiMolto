#ifndef SINGLEGROUPTRANSPORT_H
#define SINGLEGROUPTRANSPORT_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Material.h"
#include "Materials.h"
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
  int energyGroup;
  // public functions
  SingleGroupTransport(int myEnergyGroup,\
    MultiGroupTransport * myMGT,\
    Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput);

  private:
  MultiGroupTransport * MGT;
  Materials * mats;
  YAML::Node * input;
  Mesh * mesh;

};

//==============================================================================


#endif
