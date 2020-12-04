#ifndef SINGLEGROUPTRANSPORT_H
#define SINGLEGROUPTRANSPORT_H

#include "Mesh.h"
#include "Materials.h"
#include "GreyGroupQD.h"
#include "QuasidiffusionSolver.h"

using namespace std; 

class MultiGroupTransport; // forward declaration

//==============================================================================
//! SingleGroupTransport class that holds multigroup transport information

class SingleGroupTransport
{
  public:
    int energyGroup;
    arma::cube aFlux;
    arma::cube aHalfFlux;
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
    Eigen::MatrixXd calcMPQDSource(); 
    Eigen::MatrixXd calcSteadyStateMPQDSource(); 
    double calcFlux();
    double calcAlpha(string calcType="");
    double calcSteadyStateAlpha();
    void writeFlux();

  private:
    MultiGroupTransport * MGT;
    Materials * mats;
    YAML::Node * input;
    Mesh * mesh;

};

//==============================================================================


#endif
