#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "Mesh.h"
#include "Materials.h"
#include "Material.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 

//==============================================================================
//! QDSolver class that solves RZ neutron transport

class QDSolver
{
  public:
  // public functions
  QDSolver(Mesh * myMesh,\
    Materials * myMaterials,\
    YAML::Node * myInput);

  private:
  // private functions
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;
};

//==============================================================================
