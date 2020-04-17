#ifndef MultiGroupPrecursor_H
#define MultiGroupPrecursor_H

#include "Mesh.h"
#include "Materials.h"

using namespace std;

class MultiPhysicsCoupledQD;

//==============================================================================
//! Container for all precursor group

class MultiGroupPrecursor
{
  public:
  // Define default six group delayed neutron precursor data
  Eigen::VectorXd betas;
  Eigen::VectorXd lambdas;
  double beta;
  MultiGroupPrecursor(Materials * myMats,\
    Mesh * myMesh,\
    YAML::Node * myInput,\
    MultiPhysicsCoupledQD * myMPQD);
  void readInput();

  private:
  Materials * mats;
  YAML::Node * input;
  Mesh * mesh; 
  MultiPhysicsCoupledQD * mpqd;
  
};

//==============================================================================

#endif
