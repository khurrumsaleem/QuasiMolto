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
  
  // functions to map grid indices to global index
  int CFluxIndex(int iR,int iZ,int energyGroup);
  int WFluxIndex(int iR,int iZ,int energyGroup);
  int EFluxIndex(int iR,int iZ,int energyGroup);
  int NFluxIndex(int iR,int iZ,int energyGroup);
  int SFluxIndex(int iR,int iZ,int energyGroup);
  int WCurrentIndex(int iR,int iZ,int energyGroup);
  int ECurrentIndex(int iR,int iZ,int energyGroup);
  int NCurrentIndex(int iR,int iZ,int energyGroup);
  int SCurrentIndex(int iR,int iZ,int energyGroup);
  vector<int> getIndices(int iR,int iZ,int energyGroup);
  
  // functions to calculate geometry parameters
  double calcVolAvgR(double rDown,double rUp);
  vector<double> calcGeoParams(int iR,int iZ);
  
  // functions to enforce governing equations
  void assertZerothMoment(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertFirstMomentSouth(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertFirstMomentNorth(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertFirstMomentWest(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertFirstMomentEast(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertWFluxBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertEFluxBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertNFluxBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertSFluxBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertWCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertECurrentBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertNCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertSCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  double calcScatterAndFissionCoeff(int iR,int iZ,int toEnergyGroup,\
    int fromEnergyGroup);
  double calcIntegratingFactor(int iR,int iZ,double rEval,SingleGroupQD * SGQD);
  
  // function to solve linear system
  void solve();

  // function to parse solution vector
  void getFlux(SingleGroupQD * SGQD);
  Eigen::VectorXd getSolutionVector(SingleGroupQD * SGQD);
  
  // public variables
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd x;
  Eigen::VectorXd xPast;
  Eigen::VectorXd b;
  int energyGroups,nR,nZ,nGroupUnknowns;

  private:
  // private variables
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;
  // indices for accessing index and geometry parameters vectors
  int iCF = 0;
  int iWF = 1, iEF = 2, iNF = 3, iSF = 4;
  int iWC = 5, iEC = 6, iNC = 7, iSC = 8;
};

//==============================================================================
