#ifndef MultiGroupQDToMultiPhysicsQDCoupling_H
#define MultiGroupQDToMultiPhysicsQDCoupling_H

#include "Mesh.h"
#include "Materials.h"
#include "MultiPhysicsCoupledQD.h"
#include "MultiGroupQD.h"
#include "SingleGroupQD.h"
#include "HeatTransfer.h"
#include "GreyGroupQD.h"
#include "WriteData.h"

using namespace std;

//==============================================================================
//! Container for all precursor group

class MGQDToMPQDCoupling
{
  public:
    MGQDToMPQDCoupling(Mesh * myMesh,\
        Materials * myMats,\
        YAML::Node * myInput,\
        MultiPhysicsCoupledQD * myMPQD,\
        MultiGroupQD * myMGQD);
    
    // Class variables 
    double biasEps = 1E-25;
    
    // Class functions
    void initCollapsedNuclearData();
    void collapseNuclearData();
    void calculateFluxWeightedData();
    void calculateFluxWeightedInterfaceEddingtons();
    void calculateFluxWeightedBCData();
    void calculateAxialCurrentWeightedData();
    void calculateRadialCurrentWeightedData();
    void calculateRadialZetaFactors();
    void calculateAxialZetaFactors();
    void calculateCollapsedG();
    void calculateCollapsedIntFactorCoeffs();
    double checkForZeroFlux(int iZ,int iR);
    double checkForZeroInwardFluxSouthBC(int iR);
    double checkForZeroInwardFluxNorthBC(int iR);
    double checkForZeroInwardFluxEastBC(int iZ);
    double checkForZeroRadialFlux(int iZ,int iR);
    double checkForZeroAxialFlux(int iZ,int iR);
    double checkForZeroRadialCurrent(int iZ,int iR);
    double checkForZeroAxialCurrent(int iZ,int iR);
    vector<double> calcResidual(Eigen::VectorXd vector1, Eigen::VectorXd vector2);

  private:
    Mesh * mesh; 
    Materials * mats;
    YAML::Node * input;
    MultiPhysicsCoupledQD * mpqd;
    MultiGroupQD * mgqd;
};

//==============================================================================

#endif
