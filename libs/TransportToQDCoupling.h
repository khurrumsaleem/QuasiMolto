#ifndef TRANSPORTTOQDCOUPLING_H
#define TRANSPORTTOQDCOUPLING_H

#include <iostream>
#include <vector>
#include <armadillo>
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std;
using namespace arma;


class MultiGroupTransport; // forward declaration
class MultiGroupQD; // forward declaration
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
  void calcEddingtonFactors();


  private:
  YAML::Node * input;
  Mesh * mesh;
  Materials * materials;

};

//==============================================================================

#endif
