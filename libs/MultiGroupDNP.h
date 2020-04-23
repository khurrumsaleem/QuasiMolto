#ifndef MultiGroupDNP_H
#define MultiGroupDNP_H

#include "Mesh.h"
#include "Materials.h"

using namespace std;

class MultiPhysicsCoupledQD;
class SingleGroupDNP; 

//==============================================================================
//! Container for all precursor group

class MultiGroupDNP
{
  public:
  // Define default six group delayed neutron precursor data
  double beta;
  vector< shared_ptr<SingleGroupDNP> > DNPs; 
  Eigen::VectorXd recircb;
  Eigen::SparseMatrix<double> recircA;
  MultiGroupDNP(Materials * myMats,\
    Mesh * myMesh,\
    YAML::Node * myInput,\
    MultiPhysicsCoupledQD * myMPQD);
  void readInput();
  MultiPhysicsCoupledQD * mpqd;

  private:
  Materials * mats;
  YAML::Node * input;
  Mesh * mesh; 
  
};

//==============================================================================

#endif
