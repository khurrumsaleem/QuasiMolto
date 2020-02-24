#ifndef STARTINGANGLE_H
#define STARTINGANGLE_H

#include "Mesh.h"
#include "Materials.h"

using namespace std; 

//==============================================================================
//! StartingAngle class that solves RZ neutron transport at the starting angles

class StartingAngle
{
  public:
  // public functions
  StartingAngle(Mesh * myMesh,
          Materials * myMaterials,\
          YAML::Node * myInput);
  void calcStartingAngle(cube * halfAFlux,\
    Eigen::MatrixXd * source,\
    Eigen::MatrixXd * alpha,\
    int energyGroup);
  void solveAngularFlux(cube * halfAFlux,\
    Eigen::MatrixXd * source,\
    Eigen::MatrixXd * alpha,\
    int energyGroup,\
    int iXi);
  
  // default boundary conditions; homogeneous
  vector<double> upperBC,lowerBC,outerBC;
  Eigen::MatrixXd calckR(double myGamma);
  Eigen::MatrixXd calckZ(double myGamma);
  Eigen::MatrixXd calclR(double myGamma);
  Eigen::MatrixXd calclZ(double myGamma);
  Eigen::MatrixXd calct1(double myGamma);
  Eigen::MatrixXd calct2(double myGamma);
  Eigen::VectorXd calcSubCellVol(int myiZ, int myiR);
  Eigen::VectorXd calcMMSSource(int myiZ,int myiR,\
    int energyGroup,int iXi,Eigen::MatrixXd sigT, Eigen::VectorXd subCellVol);

  private:
  // private functions
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;
};

//==============================================================================

#endif
