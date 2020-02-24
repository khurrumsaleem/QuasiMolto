#ifndef MMS_H
#define MMS_H

#include "Mesh.h"
#include "MultiGroupTransport.h"
#include "MultiGroupTransport.h"
#include "SingleGroupTransport.h"
#include "SimpleCornerBalance.h"
#include "StartingAngle.h"

using namespace std; 

//==============================================================================
//! SimpleCornerBalance class that solves RZ neutron transport

class MMS
{
  public:
  // public functions
  MMS( MultiGroupTransport * myMGT,\
    Mesh * myMesh,\
    Materials * myMaterials,\
    YAML::Node * myInput);
  Eigen::MatrixXd isotropicTransportSourceMMS(double xi,\
    double mu,\
    double t);
  Eigen::MatrixXd anisotropicTransportSourceMMS(double xi,\
    double mu,\
    double t);
  void timeDependent();
  double c = 1.0;

  private:
  // private functions
  MultiGroupTransport * MGT;
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;
};

//==============================================================================

#endif
