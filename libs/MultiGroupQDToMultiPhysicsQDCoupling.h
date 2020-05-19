#ifndef MultiGroupQDToMultiPhysicsQDCoupling_H
#define MultiGroupQDToMultiPhysicsQDCoupling_H

#include "Mesh.h"
#include "Materials.h"
#include "MultiPhysicsCoupledQD.h"
#include "MultiGroupQD.h"
#include "SingleGroupQD.h"
#include "GreyGroupQD.h"

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
  void initCollapsedNuclearData();
  void solveMGQD();
  void collapseNuclearData();
  void calculateFluxWeightedData();
  void calculateFluxWeightedBCData();
  void calculateAxialCurrentWeightedData();
  void calculateRadialCurrentWeightedData();
  void calculateRadialZetaFactors();
  void calculateAxialZetaFactors();
  double checkForZeroRadialCurrent(int iZ,int iR);
  double checkForZeroAxialCurrent(int iZ,int iR);

  private:
    Mesh * mesh; 
    Materials * mats;
    YAML::Node * input;
    MultiPhysicsCoupledQD * mpqd;
    MultiGroupQD * mgqd;
};

//==============================================================================

#endif
