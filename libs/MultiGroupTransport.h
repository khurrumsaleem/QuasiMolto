#ifndef MULTIGROUPTRANSPORT_H
#define MULTIGROUPTRANSPORT_H

#include "Mesh.h"
#include "Materials.h"
#include "MultiPhysicsCoupledQD.h"
#include "MultiGroupQD.h"

using namespace std; 

class SingleGroupTransport; // forward declaration
class StartingAngle; // forward declaration
class SimpleCornerBalance; // forward declaration

//==============================================================================
//! MultiGroupTransport class that holds multigroup transport information

class MultiGroupTransport
{
  public:
    vector< shared_ptr<SingleGroupTransport> > SGTs;
    shared_ptr<StartingAngle> startAngleSolve;
    shared_ptr<SimpleCornerBalance> SCBSolve;
    // public functions
    MultiGroupTransport(Materials * myMaterials,\
        Mesh * myMesh,\
        YAML::Node * myInput);

    // default convergence criteria
    double epsAlpha = 1E-3;
    double epsFlux = 1E-5;
    double epsFissionSource = 1E-5;

    // defaults for maximum iterations
    double sourceMaxIter = 500;
    double powerMaxIter = 10000;

    // Boolean to determine use of grey group sources
    bool useMPQDSources = false;    
 
    // Pointers
    MultiPhysicsCoupledQD * mpqd;
    MultiGroupQD * mgqd;

    // public functions
    void solveStartAngles();
    void solveSCBs();
    bool calcSources(string calcType = "FS");
    bool calcFluxes(string printResidual = "noprint");
    bool calcAlphas(string printResidual = "noprint", string calcType = "");
    bool calcFissionSources(string printResidual = "noprint");
    bool sourceIteration();
    bool powerIteration();
    void solveTransportOnly();
    void printDividers();
    void writeFluxes();
    void assignMultiPhysicsCoupledQDPointers(MultiGroupQD * myMGQD,\
        MultiPhysicsCoupledQD * myMPQD);

  private:
    YAML::Node * input;
    Mesh * mesh;
    Materials * materials;
    // for print formatting
    int spacing = 15;

};

//==============================================================================


#endif
