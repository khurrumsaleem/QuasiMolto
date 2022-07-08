#ifndef TRANSPORTTOQDCOUPLING_H
#define TRANSPORTTOQDCOUPLING_H

#include "Mesh.h"
#include "MultiGroupQD.h"
#include "MultiGroupTransport.h"
#include "SingleGroupQD.h"
#include "SingleGroupTransport.h"
#include "QuasidiffusionSolver.h"

using namespace std;

//==============================================================================
//! TransportToQDCoupling class that handles communication between transport and 
///   quasidiffusion objects   

class TransportToQDCoupling
{
  public:
  MultiGroupTransport * MGT;
  MultiGroupQD * MGQD;
  TransportToQDCoupling(Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput,\
    MultiGroupTransport * myMGT,\
    MultiGroupQD * myMGQD);
 
  // User-definable variables
  double epsEddington = 1.0E-5;

  // Class functions
  bool calcEddingtonFactors();
  bool calcInterfaceEddingtonFactors();
  void calcGFactors();
  void calcIntFactorCoeffs();
  void calcBCs();
  void updateTransportFluxes();
  void updateTransportPrevFluxes();
  double calcResidual(Eigen::MatrixXd matrix1,Eigen::MatrixXd matrix2);

  // Utility functions
  void checkOptionalParams();

  private:
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;

};

//==============================================================================

#endif
