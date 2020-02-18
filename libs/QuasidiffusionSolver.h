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
  vector<int> getIndices(int iR,int iZ,int energyGroup);
  
  // functions to calculate geometry parameters
  double calcVolAvgR(double rDown,double rUp);
  vector<double> calcGeoParams(int iR,int iZ);
  
  // functions to enforce governing equations
  void assertZerothMoment(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void applyRadialBoundary(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void applyAxialBoundary(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  
  // functions to enforce coefficients for facial currents
  void southCurrent(double coeff,int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void northCurrent(double coeff,int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void westCurrent(double coeff,int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void eastCurrent(double coeff,int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);

  // functions to assert flux BCs
  void assertWFluxBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertEFluxBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertNFluxBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertSFluxBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);

  // functions to assert current BCs
  void assertWCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertECurrentBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertNCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertSCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);

  // wrapper to assert either a flux or current BC
  void assertWBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertEBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertNBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);
  void assertSBC(int iR,int iZ,int iEq,int energyGroup,\
    SingleGroupQD * SGQD);

  double calcScatterAndFissionCoeff(int iR,int iZ,int toEnergyGroup,\
    int fromEnergyGroup);
  double calcIntegratingFactor(int iR,int iZ,double rEval,SingleGroupQD * SGQD);
  
  // function to solve linear system
  void solve();

  // function to parse solution vector
  void getFlux(SingleGroupQD * SGQD);
  Eigen::VectorXd getFluxSolutionVector(SingleGroupQD * SGQD);
  Eigen::VectorXd getCurrentSolutionVector(SingleGroupQD * SGQD);

  // check for input parameters
  void checkOptionalParams();
  
  // public variables
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd x;
  Eigen::VectorXd xPast,currPast;
  Eigen::VectorXd b;
  int energyGroups,nR,nZ,nGroupUnknowns,nGroupCurrentUnknowns;
  bool reflectingBCs = false;
  
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
