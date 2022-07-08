#ifndef MultiPhysicsCoupledQD_H
#define MultiPhysicsCoupledQD_H

#include "Mesh.h"
#include "Materials.h"
#include "GreyGroupSolver.h"
#include "SingleGroupDNP.h"
#include "WriteData.h"

using namespace std;

class HeatTransfer;
class MultiGroupDNP;
class GreyGroupQD;

//==============================================================================
//! Contains precursor, heat, and grey group quasidiffusion objects.
//!   Also owns and is responsible for solving the coupled system each of
//!  of those object contribute to.

class MultiPhysicsCoupledQD
{
  public:
    // Constructor
    MultiPhysicsCoupledQD(Materials * myMats,\
        Mesh * myMesh,\
        YAML::Node * myInput);

    // Pointers
    HeatTransfer * heat;
    MultiGroupDNP * mgdnp;
    GreyGroupQD * ggqd;

    // Class variables
    int nUnknowns;
    string outputDir = "MPQD/";
    
    // User-definable variables
    double epsFlux = 1E-6, epsTemp = 1E-6;
    
    // MPQD linear system variables
    Vec x_p,xPast_p,b_p;
    Vec xPast_p_seq;
    Mat A_p;
    KSP ksp;
    PC pc;

    // Dual purpose
    int solve();

    // Steady state 
    int buildSteadyStateLinearSystem();
    void updateSteadyStateVarsAfterConvergence();
    void solveSteadyState();

    // Transient
    int buildLinearSystem();
    void updateVarsAfterConvergence();
    void solveTransient();

    // Pseudo transient
    int buildPseudoTransientLinearSystem();
    void updatePseudoTransientVars();

    // Utility functions 
    int  fluxSource(int iZ, int iR, int iEq, double coeff);
    int  dnpSource(int iZ, int iR, int iEq, double coeff);
    void setInitialCondition();
    void writeVars();
    void printVars();
    void checkOptionalParams();

  private:
    Materials * mats;
    Mesh * mesh;
    YAML::Node * input;

};

//==============================================================================

#endif
