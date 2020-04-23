#ifndef MultiPhysicsCoupledQD_H
#define MultiPhysicsCoupledQD_H

#include "Mesh.h"
#include "Materials.h"

using namespace std;

class HeatTransfer;
class MultiGroupDNP;

//==============================================================================
//! Contains precursor, heat, and grey group quasidiffusion objects.
//!   Also owns and is responsible for solving the coupled system each of
//!  of those object contribute to.

class MultiPhysicsCoupledQD
{
  public:
  MultiPhysicsCoupledQD(Materials * myMats,\
    Mesh * myMesh,\
    YAML::Node * myInput);
  
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd x,b;
  void fluxSource(int iZ,int iR,int iEq,double coeff);
  void solveLinearSystem();
  HeatTransfer * heat;
  MultiGroupDNP * mgdnp;

  private:
  Materials * mats;
  Mesh * mesh;
  YAML::Node * input;

};

//==============================================================================

#endif
