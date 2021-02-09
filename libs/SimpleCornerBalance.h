#ifndef SIMPLECORNERBALANCE_H
#define SIMPLECORNERBALANCE_H

#include "Mesh.h"
#include "Materials.h"

using namespace std; 

//==============================================================================
//! SimpleCornerBalance class that solves RZ neutron transport

class SimpleCornerBalance
{
  public:
  // public functions
  SimpleCornerBalance(Mesh * myMesh,\
    Materials * myMaterials,\
    YAML::Node * myInput);
  void solve(arma::cube * aFlux,\
    arma::cube * halfAFlux,\
    Eigen::MatrixXd * source,\
    Eigen::MatrixXd * alpha,\
    int energyGroup);
  void solveAngularFlux(arma::cube * aFlux,\
    arma::cube * halfAFlux,\
    Eigen::MatrixXd * source,\
    Eigen::MatrixXd * alpha,\
    int energyGroup,int iXi,int iMu);
  void solveAngularFluxNegMu(arma::cube * aFlux,\
    arma::cube * halfAFlux,\
    Eigen::MatrixXd * source,\
    Eigen::MatrixXd * alpha,\
    int energyGroup,int iXi,int iMu);
  void solveAngularFluxPosMu(arma::cube * aFlux,\
    arma::cube * halfAFlux,\
    Eigen::MatrixXd * source,\
    Eigen::MatrixXd * alpha,\
    int energyGroup,int iXi,int iMu);
  
  // default boundary conditions; homogeneous
  vector<double> upperBC,lowerBC,outerBC;
  Eigen::MatrixXd calckR(double myGamma);
  Eigen::MatrixXd calckZ(double myGamma);
  Eigen::MatrixXd calclR(double myGamma);
  Eigen::MatrixXd calclZ(double myGamma);
  Eigen::MatrixXd calct(double myGamma);
  Eigen::MatrixXd calcR(double myGamma);
  Eigen::VectorXd calcSubCellVol(int myiZ, int myiR);
  Eigen::VectorXd calcMMSSource(int myiZ,int myiR,int energyGroup,\
    int iXi, int iMu, Eigen::MatrixXd sigT, Eigen::VectorXd subCellVol);


  private:
  // private functions
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;
};

//==============================================================================

#endif
