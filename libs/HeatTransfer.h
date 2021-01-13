#ifndef HeatTransfer_H
#define HeatTransfer_H 

#include "Mesh.h"
#include "Materials.h"
#include "GreyGroupQD.h"

using namespace std;

class MultiPhysicsCoupledQD;

//==============================================================================
//! Contains information and builds linear system for heat transfer

class HeatTransfer
{
  public:
  HeatTransfer(Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput,\
    MultiPhysicsCoupledQD * myQD);

  // Default temperatures taken from "Introduction to Moltres:..." [2018]
  double wallT = 922.0;
  double inletT = 922.0;
  int coreInletIndex,coreOutletIndex,nUnknowns,indexOffset = 0; 
  string modIrradiation = "volume", axial = "axial", volume = "volume",\
                           fuel = "fuel";
  string fluxLimiter = "superbee";
  Eigen::SparseMatrix<double,Eigen::RowMajor> Atemp;
  Eigen::MatrixXd temp,flux,dirac,inletTemp;
  Eigen::VectorXd inletDensity,inletVelocity,inletcP,outletTemp;        
  int getIndex(int iZ,int iR);
  void buildLinearSystem();
  void buildSteadyStateLinearSystem();
  void gammaSource(int iZ,int iR,int iEq,double coeff,\
      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> * myA);
  Eigen::MatrixXd calcExplicitFissionEnergy();
  Eigen::MatrixXd calcExplicitAxialFissionEnergy();
  Eigen::MatrixXd calcExplicitAxialFuelFissionEnergy();
  void calcDiracs();
  void calcFluxes();
  void calcImplicitFluxes();
  void getTemp();
  void setTemp();
  Eigen::MatrixXd returnCurrentTemp();
  void assignBoundaryIndices();
  void updateBoundaryConditions();
  double calcPhi(double theta,string fluxLimiter);
  double calcTheta(double TupwindInterface,double Tinterface);
  void checkOptionalParams();

  /* PETSc functions */

  // Steady state functions

  int buildSteadyStateLinearSystem_p();
  
  // Transient functions
  
  private:
  Materials * mats;
  Mesh * mesh;
  YAML::Node * input;
  MultiPhysicsCoupledQD * mpqd;
  // Indices for volume, west face surface area, east face surface area, etc...
  int iVol = 0;
  int iWF = 1, iEF = 2, iNF = 3, iSF = 4;
};

//==============================================================================

#endif
