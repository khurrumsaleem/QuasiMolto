#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "Mesh.h"
#include "Materials.h"
#include "Material.h"
#include "SingleGroupQD.h"
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
  void formLinearSystem(SingleGroupQD * SGQD);
  int CFluxIndex(int iR,int iZ,int energyGroup);
  int WFluxIndex(int iR,int iZ,int energyGroup);
  int EFluxIndex(int iR,int iZ,int energyGroup);
  int NFluxIndex(int iR,int iZ,int energyGroup);
  int SFluxIndex(int iR,int iZ,int energyGroup);
  int WCurrentIndex(int iR,int iZ,int energyGroup);
  int ECurrentIndex(int iR,int iZ,int energyGroup);
  int NCurrentIndex(int iR,int iZ,int energyGroup);
  int SCurrentIndex(int iR,int iZ,int energyGroup);
  double calcVolAvgR(double rDown,double rUp);
  vector<int> indices(int iR,int iZ,int energyGroup);
  vector<double> calcGeoParams(int iR,int iZ);
  
  // public variables
  Eigen::SparseMatrix<double> A;
  int energyGroups,nR,nZ,nGroupUnknowns;

  private:
  // private variables
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;
  // indices for accessing index vector
  int iCF = 0;
  int iWF = 1, iEF = 2, iNF = 3, iSF = 4;
  int iWC = 5, iEC = 6, iNC = 7, iSC = 8;
};

//==============================================================================
