#ifndef SingleGroupDNP_H
#define SingleGroupDNP_H

#include "Mesh.h"

using namespace std;

class MultiGroupDNP; // forward declaration

//==============================================================================
//! Contains information and builds linear system for a single precursor group

class SingleGroupDNP
{
  public:
  SingleGroupDNP(MultiGroupDNP * myMGDNPS,\
    double myBeta,\
    double myLambda);
  double beta; 
  double lambda; 

  private: 
  MultiGroupDNP * mgdnp;
};

//==============================================================================

#endif
