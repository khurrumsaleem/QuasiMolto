#include "Mesh.h"
#include "Materials.h"
#include "Material.h"
#include "MultiGroupTransport.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"

using namespace std; 

//==============================================================================
//! SimpleCornerBalance class that solves RZ neutron transport

class MMS
{
  public:
  // public functions
  MMS( MultiGroupTransport * myMGT,\
    Mesh * myMesh,\
    Materials * myMaterials,\
    YAML::Node * myInput);
  Eigen::MatrixXd isotropicTransportSourceMMS(double xi,\
    double mu,\
    double t);
  void timeDependent();
  double c = 1.0;

  private:
  // private functions
  MultiGroupTransport * MGT;
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;
};

//==============================================================================
