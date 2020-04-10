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
  
  private:
  Materials * mats;
  Mesh * mesh;
  YAML::Node * input;
  MultiPhysicsCoupledQD * qd;
};

//==============================================================================

#endif
