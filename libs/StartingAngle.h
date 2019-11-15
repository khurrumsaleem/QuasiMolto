#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "Mesh.h"
#include "Materials.h"
#include "Material.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

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
  
  // default boundary conditions; homogeneous
  double upperBC=0,lowerBC=0,outerBC=0;
  Eigen::MatrixXd calckR(double myGamma);
  Eigen::MatrixXd calckZ(double myGamma);
  Eigen::MatrixXd calclR(double myGamma);
  Eigen::MatrixXd calclZ(double myGamma);
  Eigen::MatrixXd calct1(double myGamma);
  Eigen::MatrixXd calct2(double myGamma);
  Eigen::VectorXd calcSubCellVol(int myiZ, int myiR);

  private:
  // private functions
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;
};

//==============================================================================
