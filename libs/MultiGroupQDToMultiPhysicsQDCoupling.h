#ifndef MultiGroupQDToMultiPhysicsQDCoupling_H
#define MultiGroupQDToMultiPhysicsQDCoupling_H

#include "Mesh.h"
#include "Materials.h"
#include "MultiPhysicsCoupledQD.h"
#include "MultiGroupQD.h"

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

  private:
    Mesh * mesh; 
    Materials * mats;
    YAML::Node * input;
    MultiPhysicsCoupledQD * mpqd;
    MultiGroupQD * mgqd;
};

//==============================================================================

#endif
