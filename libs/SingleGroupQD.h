#ifndef SINGLEGROUPQD_H 
#define SINGLEGROUPQD_H

#include <iostream>
#include <fstream>
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

class MultiGroupQD; // forward declaration

//==============================================================================
//! SingleGroupQD class that holds quasidiffusion information

class SingleGroupQD
{
  public:
  int energyGroup;
  Eigen::MatrixXd sFlux;
  Eigen::MatrixXd sFluxPrev;
  Eigen::MatrixXd q;
  Eigen::MatrixXd fissionSource;
  Eigen::MatrixXd scatterSource;
  // public functions
  SingleGroupQD(int myEnergyGroup,\
    MultiGroupQD * myMGQD,\
    Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput);

  private:
  MultiGroupQD * MGQD;
  Materials * mats;
  YAML::Node * input;
  Mesh * mesh;

};

//==============================================================================

#endif
