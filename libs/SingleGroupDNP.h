#ifndef SingleGroupDNP_H
#define SingleGroupDNP_H

#include "Mesh.h"
#include "Materials.h"

using namespace std;

class MultiGroupDNP; // forward declaration

//==============================================================================
//! Contains information and builds linear system for a single precursor group

class SingleGroupDNP
{
  public:
  Eigen::MatrixXd dnpConc,recircConc,flux,recircFlux,dirac,recircDirac;
  Eigen::MatrixXd inletConc,recircInletConc;
  Eigen::VectorXd inletVelocity,recircInletVelocity,outletConc,recircOutletConc;
  double beta; 
  double lambda;
  int coreInletIndex,coreOutletIndex,recircInletIndex,recircOutletIndex;
  string fluxLimiter = "superbee";
  int coreIndexOffset = 0;
  int recircIndexOffset = 0;
  SingleGroupDNP(Materials * myMats,\
    Mesh * myMesh,\
    MultiGroupDNP * myMGDNPS,\
    double myBeta,\
    double myLambda);
  int getIndex(int iZ,int iR,int indexOffset);
  void assignBoundaryIndices();
  void updateBoundaryConditions();
  void calcRecircDNPFluxes();
  void calcCoreDNPFluxes();
  void buildRecircLinearSystem();
  void buildCoreLinearSystem();
  void buildLinearSystem(Eigen::SparseMatrix<double> * myA,\
    Eigen::VectorXd * myb,\
    Eigen::MatrixXd myDNPConc,\
    Eigen::MatrixXd myDNPFlux,\
    Eigen::MatrixXd mySigF,\
    rowvec dzs,\
    int myIndexOffset,
    bool fluxSource = true);
  Eigen::MatrixXd calcDiracs(Eigen::MatrixXd dnpConc,\
    Eigen::MatrixXd inletConc,\
    Eigen::VectorXd outletConc);
  Eigen::MatrixXd calcFluxes(Eigen::MatrixXd myDNPConc,\
    Eigen::MatrixXd myFlowVelocity,\
    Eigen::MatrixXd dirac,\
    Eigen::MatrixXd inletConc,\
    Eigen::VectorXd inletVelocity,\
    rowvec dzs);
  void getCoreConc();
  void getRecircConc();
  double calcPhi(double theta,string fluxLimiter); 
  double calcTheta(double DNPupwindInterface,double DNPinterface);

  private: 
  Materials * mats;
  Mesh * mesh;
  MultiGroupDNP * mgdnp;
};

//==============================================================================

#endif
