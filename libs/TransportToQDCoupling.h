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
  bool calcEddingtonFactors();
  void calcBCs();
  void solveTransportWithQDAcceleration();
  void updateTransportFluxes();
  void updateTransportPrevFluxes();
  void checkOptionalParams();
  double epsEddington = 1.0E-5;


  private:
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;

};

//==============================================================================

#endif
