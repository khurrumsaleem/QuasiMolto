#ifndef SINGLEGROUPTRANSPORT_H
#define SINGLEGROUPTRANSPORT_H

#include "Mesh.h"
#include "Materials.h"

using namespace std; 
using namespace arma;

class MultiGroupTransport; // forward declaration

//==============================================================================
//! SingleGroupTransport class that holds multigroup transport information

class SingleGroupTransport
{
  public:
  int energyGroup;
  cube aFlux;
  cube aHalfFlux;
  Eigen::MatrixXd sFlux; 
  Eigen::MatrixXd sFluxPrev; 
  Eigen::MatrixXd alpha; 
  Eigen::MatrixXd q; 
  Eigen::MatrixXd fissionSource; 
  Eigen::MatrixXd scatterSource; 
  // public functions
  SingleGroupTransport(int myEnergyGroup,\
    MultiGroupTransport * myMGT,\
    Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput);
  void solveStartAngle();
  void solveSCB();
  double calcSource(string calcType="FS"); 
  double calcFissionSource(); 
  double calcScatterSource(); 
  double calcFlux();
  double calcAlpha();
  void writeFlux();

  private:
  MultiGroupTransport * MGT;
  Materials * mats;
  YAML::Node * input;
  Mesh * mesh;

};

//==============================================================================


#endif
