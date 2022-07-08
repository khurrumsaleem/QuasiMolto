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
#include "PETScWrapper.h"

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
    
    // Class variables 
    string outputDir = "Solve_Metrics/";
   
    // User-definable variables
    double resetThreshold = 1E100, relaxTolELOT = 3E-4, relaxTolMGLOQD = 3E-4,\
           ratedPower = 8e6, fluxNormalization = 1, epsK = 1E-8;
    bool skipMGHOT = false;
    
    // Steady state
    void solveSteadyStateResidualBalance(bool outputVars);
    void solveSteadyStateMGHOT();
    void solveSteadyStateMGLOQD();
    void solveSteadyStateELOT();

    // Transient
    bool solveOneStepResidualBalance(bool outputVars);
    void solveTransient();
    void solveSteadyStateTransientResidualBalance(bool outputVars);
    void solveMGHOT();
    void solveMGLOQD();
    void solveELOT();

    // Pseudo transient
    void solveSteadyStatePseudoTransient(bool outputVars);
    void solvePseudoTransient();
    bool solvePseudoTransientResidualBalance(bool outputVars);
    void solvePseudoTransientELOT();

    // Utility functions
    double calcK(Eigen::MatrixXd oldFlux, Eigen::MatrixXd newFlux,\
        Eigen::MatrixXd volume, double kold);
    double eps(double residual, double relaxationTolerance = 1E-14);
    double epsTemp(double residual, double relaxationTolerance = 1E-14);
    double relaxedEpsK(double residual, double relaxationTolerance = 1E-14);
    void   checkOptionalParameters();

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
