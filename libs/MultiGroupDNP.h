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
  int indexOffset = 0;
  int nCoreUnknowns,nRecircUnknowns;
  double beta;
  vector< shared_ptr<SingleGroupDNP> > DNPs; 
  Eigen::VectorXd recircb,recircx;
  Eigen::SparseMatrix<double> recircA;
  MultiGroupDNP(Materials * myMats,\
    Mesh * myMesh,\
    YAML::Node * myInput,\
    MultiPhysicsCoupledQD * myMPQD,\
    int myIndexOffset);
  void readInput();
  void buildRecircLinearSystem();
  void buildCoreLinearSystem();
  void getCoreDNPConc();
  void setCoreDNPConc();
  void printCoreDNPConc();
  void solveRecircLinearSystem();
  MultiPhysicsCoupledQD * mpqd;
  YAML::Node * input;

  private:
  Materials * mats;
  Mesh * mesh; 
  
};

//==============================================================================

#endif
