#ifndef SingleGroupPrecursor_H
#define SingleGroupPrecursor_H

#include "Mesh.h"

using namespace std;

class MultiGroupPrecursor; // forward declaration

//==============================================================================
//! Contains information and builds linear system for a single precursor group

class SingleGroupPrecursor
{
  public:
  SingleGroupPrecursor(MultiGroupPrecursor * myMGP,\
    double myBeta,\
    double myLambda);
  double beta; 
  double lambda; 

  private: 
  MultiGroupPrecursor * dnps;
};

//==============================================================================

#endif
