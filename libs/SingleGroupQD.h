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
  Eigen::MatrixXd sFluxR;
  Eigen::MatrixXd sFluxZ;
  Eigen::MatrixXd currentR;
  Eigen::MatrixXd currentZ;
  Eigen::MatrixXd sFluxPrev;
  Eigen::MatrixXd q;
  Eigen::MatrixXd fissionSource;
  Eigen::MatrixXd scatterSource;
  
  // Eddington factors
  Eigen::MatrixXd Err;
  Eigen::MatrixXd Ezz;
  Eigen::MatrixXd Erz;

  // flux boundary conditions  
  Eigen::VectorXd wFluxBC;
  Eigen::VectorXd eFluxBC;
  Eigen::VectorXd nFluxBC;
  Eigen::VectorXd sFluxBC;
  
  // current boundary conditions
  Eigen::VectorXd wCurrentRBC;
  Eigen::VectorXd eCurrentRBC;
  Eigen::VectorXd nCurrentZBC;
  Eigen::VectorXd sCurrentZBC;

  // public functions
  SingleGroupQD(int myEnergyGroup,\
    MultiGroupQD * myMGQD,\
    Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput);
  void formContributionToLinearSystem();
  void getFlux();
  Eigen::VectorXd getSolutionVector();

  private:
  MultiGroupQD * MGQD;
  Materials * mats;
  YAML::Node * input;
  Mesh * mesh;

};

//==============================================================================

#endif
