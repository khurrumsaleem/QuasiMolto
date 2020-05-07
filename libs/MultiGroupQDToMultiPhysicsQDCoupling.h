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
    MGQDToMPQDCoupling();

  private:
    Mesh * mesh; 
    Materials * mats;
    YAML::Node * input;
    MultiPhysicsCoupledQD * mpqd;
    MultiGroupQD * MGQD;
};

//==============================================================================

#endif
