#ifndef MultiGroupPrecursor_H
#define MultiGroupPrecursor_H

#include "Mesh.h"
#include "Materials.h"

using namespace std;

class MultiPhysicsCoupledQD;
class SingleGroupPrecursor; 

//==============================================================================
//! Container for all precursor group

class MultiGroupPrecursor
{
  public:
  // Define default six group delayed neutron precursor data
  //Eigen::VectorXd betas;
  //Eigen::VectorXd lambdas;
  double beta;
  vector< shared_ptr<SingleGroupPrecursor> > DNPs; 
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
