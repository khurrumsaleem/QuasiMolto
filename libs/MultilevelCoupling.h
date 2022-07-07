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
    double resetThreshold = 1E100, relaxTolELOT = 3E-4, relaxTolMGLOQD = 3E-4,\
           ratedPower = 8e6, fluxNormalization = 1, epsK = 1E-8;
    bool p1Approx = false, iterativeMGLOQD = false, iterativeELOT = false;
    void solveMGHOT();
    void solveSteadyStateMGHOT();
    double calcK(Eigen::MatrixXd oldFlux, Eigen::MatrixXd newFlux,\
        Eigen::MatrixXd volume, double kold);
    double eps(double residual, double relaxationTolerance = 1E-14);
    double epsTemp(double residual, double relaxationTolerance = 1E-14);
    double relaxedEpsK(double residual, double relaxationTolerance = 1E-14);
    void checkOptionalParameters();
    string outputDir = "Solve_Metrics/";
    
    /* PETSc FUNCTIONS */

    // Steady state
    void solveSteadyStateResidualBalance(bool outputVars);
    void solveSteadyStateMGLOQD();
    void solveSteadyStateELOT();

    // Transient
    bool solveOneStepResidualBalance(bool outputVars);
    void solveTransient();
    void solveSteadyStateTransientResidualBalance(bool outputVars);
    void solveMGLOQD();
    void solveELOT();

    // Pseudo transient
    void solveSteadyStatePseudoTransient(bool outputVars);
    void solvePseudoTransient();
    bool solvePseudoTransientResidualBalance(bool outputVars);
    void solvePseudoTransientELOT();


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
