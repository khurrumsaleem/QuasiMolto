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
  SingleGroupDNP(Materials * myMats,\
    Mesh * myMesh,\
    MultiGroupDNP * myMGDNPS,\
    double myBeta,\
    double myLambda);
  double beta; 
  double lambda; 

  private: 
  Materials * mats;
  Mesh * mesh;
  MultiGroupDNP * mgdnp;
};

//==============================================================================

#endif
