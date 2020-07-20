#ifndef MULTILEVELCOUPLING_H
#define MULTILEVELCOUPLING_H

#include "Mesh.h"
#include "Materials.h"
#include "MultiPhysicsCoupledQD.h"
#include "MultiGroupQD.h"
#include "SingleGroupQD.h"
#include "SingleGroupTransport.h"
#include "MultiGroupTransport.h"
#include "TransportToQDCoupling.h"
#include "MultiGroupQDToMultiPhysicsQDCoupling.h"
#include "GreyGroupQD.h"
#include "WriteData.h"

using namespace std;

//==============================================================================
//!  Class to manipulate MGT, MGQD, and MPQD objects for a coupled solve

class MultilevelCoupling
{
  public:
    MultilevelCoupling(Mesh * myMesh,\
        Materials * myMats,\
        YAML::Node * myInput,\
        MultiGroupTransport * myMGT,\
        MultiGroupQD * myMGQD,\
        MultiPhysicsCoupledQD * myMPQD);
    double resetThreshold = 1E100, relaxTolELOT = 3E-4, relaxTolMGLOQD = 3E-4;
    bool solveOneStep();
    bool solveOneStepResidualBalance();
    bool initialSolve();
    void solveMGHOT();
    void solveMGLOQD();
    void solveELOT();
    void solveTransient();
    double eps(double residual, double relaxationTolerance = 1E-14);
    void checkOptionalParameters();
    string outputDir = "Solve_Metrics/";
    

  private:
    Mesh * mesh; 
    Materials * mats;
    YAML::Node * input;
    MultiGroupTransport * mgt;
    MultiPhysicsCoupledQD * mpqd;
    MultiGroupQD * mgqd;
    TransportToQDCoupling * MGTToMGQD; 
    MGQDToMPQDCoupling * MGQDToMPQD; 
};

//==============================================================================

#endif
