#ifndef HeatTransfer_H
#define HeatTransfer_H 

#include "Mesh.h"
#include "Materials.h"
#include "MultiPhysicsCoupledQD.h"

using namespace std;

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
  int coreInletIndex,coreOutletIndex,indexOffset = 0; 
  string fluxLimiter = "superbee";
  Eigen::MatrixXd temp,flux,dirac,inletTemp;
  Eigen::VectorXd inletDensity,inletVelocity,inletcP,outletTemp;        
  int getIndex(int iZ,int iR);
  void calcDiracs();
  void calcFluxes();
  void getTemp();
  void assignBoundaryIndices();
  void updateBoundaryConditions();
  double calcPhi(double theta,string fluxLimiter);
  
  private:
  Materials * mats;
  Mesh * mesh;
  YAML::Node * input;
  MultiPhysicsCoupledQD * qd;
  // Indices for volume, west face surface area, east face surface area, etc...
  int iVol = 0;
  int iWF = 1, iEF = 2, iNF = 3, iSF = 4;
};

//==============================================================================

#endif
