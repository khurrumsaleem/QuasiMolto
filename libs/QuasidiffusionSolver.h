#ifndef QUASIDIFFUSIONSOLVER_H
#define QUASIDIFFUSIONSOLVER_H

#include "Mesh.h"
#include "Materials.h"
#include "MultiPhysicsCoupledQD.h"
#include "PETScWrapper.h"

class SingleGroupQD;
class GreyGroupQD;

using namespace std; 

//==============================================================================
//! QDSolver class that solves RZ neutron transport

class QDSolver
{
  public:
    // public functions
    QDSolver(Mesh * myMesh,\
        Materials * myMaterials,\
        YAML::Node * myInput);

    // functions to map grid indices to global index
    vector<int> getIndices(int iR,int iZ,int energyGroup);

    // functions to calculate geometry parameters
    double calcVolAvgR(double rDown,double rUp);
    vector<double> calcGeoParams(int iR,int iZ);
    
    // functions to get Eddington factors
    double getWestErr(int iZ,int iR, SingleGroupQD * SGQD);
    double getWestErz(int iZ,int iR, SingleGroupQD * SGQD);
    double getEastErr(int iZ,int iR, SingleGroupQD * SGQD);
    double getEastErz(int iZ,int iR, SingleGroupQD * SGQD);
    double getNorthEzz(int iZ,int iR, SingleGroupQD * SGQD);
    double getNorthErz(int iZ,int iR, SingleGroupQD * SGQD);
    double getSouthEzz(int iZ,int iR, SingleGroupQD * SGQD);
    double getSouthErz(int iZ,int iR, SingleGroupQD * SGQD);

    double calcScatterAndFissionCoeff(int iR,int iZ,int toEnergyGroup,\
        int fromEnergyGroup);
    double calcIntegratingFactor(int iR,int iZ,double rEval,\
        SingleGroupQD * SGQD);

    // function to parse solution vector
    int getFlux(SingleGroupQD * SGQD);
    Eigen::VectorXd getFluxSolutionVector(SingleGroupQD * SGQD);
    Eigen::VectorXd getCurrentSolutionVector(SingleGroupQD * SGQD);

    // check for input parameters
    void checkOptionalParams();

    // Class variables
    int energyGroups,nR,nZ,nGroupUnknowns,nGroupCurrentUnknowns;
    int nUnknowns,nCurrentUnknowns;
    bool useMPQDSources = false;
    MultiPhysicsCoupledQD * mpqd;

    // User-definable variables
    bool reflectingBCs = false;
    bool goldinBCs = false;
    bool p1BCs = false;
    bool fluxBCs = false;

    /* PETSc stuff */
    // PETSc variables (the p denotes a petsc variable)
    Vec x_p,b_p,d_p;
    Vec xPast_p,currPast_p;
    Vec xPast_p_seq,currPast_p_seq;
    Mat A_p,C_p;
    KSP ksp;
    PC pc;
   
    /* =========== TRANSIENT AND STEADY-STATE FUNCTIONS =============*/ 
    int solve();
    int backCalculateCurrent();

    // functions to assert flux BCs
    int assertWFluxBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int assertEFluxBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int assertNFluxBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int assertSFluxBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);

    /*=================== TRANSIENT FUNCTIONS =======================*/
    
    void formLinearSystem(SingleGroupQD * SGQD);
    void formBackCalcSystem(SingleGroupQD * SGQD);

    // functions to enforce steady state governing equations
    int assertZerothMoment(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void applyRadialBoundary(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void applyAxialBoundary(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);

    // functions to enforce coefficients for steady state facial currents
    int southCurrent(double coeff,int iR,int iZ,int iEq,\
        int energyGroup, SingleGroupQD * SGQD);
    int northCurrent(double coeff,int iR,int iZ,int iEq,\
        int energyGroup, SingleGroupQD * SGQD);
    int westCurrent(double coeff,int iR,int iZ,int iEq,\
        int energyGroup, SingleGroupQD * SGQD);
    int eastCurrent(double coeff,int iR,int iZ,int iEq,\
        int energyGroup, SingleGroupQD * SGQD);

    // functions to assert steady state current BCs
    void assertWCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertECurrentBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertNCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertSCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);

    // functions to assert steady state Gol'din BCs
    int assertNGoldinBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int assertSGoldinBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int assertEGoldinBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);

    // functions to assert Gol'din diffusion BCs
    int assertNGoldinP1BC(int iR,int iZ,int iEq,\
        int energyGroup,SingleGroupQD * SGQD);
    int assertSGoldinP1BC(int iR,int iZ,int iEq,\
        int energyGroup,SingleGroupQD * SGQD);
    int assertEGoldinP1BC(int iR,int iZ,int iEq,\
        int energyGroup,SingleGroupQD * SGQD);

    // wrapper to assert either a flux or current BC
    void assertWBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertEBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertNBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertSBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);

    // functions to enforce coefficients for calculation of steady state 
    // facial currents
    int calcSouthCurrent(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int calcNorthCurrent(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int calcWestCurrent(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int calcEastCurrent(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);

    int greyGroupSources(int iR,int iZ,int iEq,int toEnergyGroup,\
        vector<double> geoParams);
    
    /*=================== STEADY-STATE FUNCTIONS ====================*/

    void formSteadyStateLinearSystem(SingleGroupQD * SGQD);
    void formSteadyStateBackCalcSystem(SingleGroupQD * SGQD);

    // functions to enforce steady state governing equations
    int assertSteadyStateZerothMoment(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void applySteadyStateRadialBoundary(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void applySteadyStateAxialBoundary(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);

    // functions to enforce coefficients for steady state facial currents
    int steadyStateSouthCurrent(double coeff,int iR,int iZ,int iEq,\
        int energyGroup, SingleGroupQD * SGQD);
    int steadyStateNorthCurrent(double coeff,int iR,int iZ,int iEq,\
        int energyGroup, SingleGroupQD * SGQD);
    int steadyStateWestCurrent(double coeff,int iR,int iZ,int iEq,\
        int energyGroup, SingleGroupQD * SGQD);
    int steadyStateEastCurrent(double coeff,int iR,int iZ,int iEq,\
        int energyGroup, SingleGroupQD * SGQD);

    // functions to assert steady state current BCs
    void assertSteadyStateWCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertSteadyStateECurrentBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertSteadyStateNCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertSteadyStateSCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);

    // functions to assert steady state Gol'din BCs
    int assertSteadyStateNGoldinBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int assertSteadyStateSGoldinBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int assertSteadyStateEGoldinBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);

    // functions to assert Gol'din diffusion BCs
    int assertSteadyStateNGoldinP1BC(int iR,int iZ,int iEq,\
        int energyGroup,SingleGroupQD * SGQD);
    int assertSteadyStateSGoldinP1BC(int iR,int iZ,int iEq,\
        int energyGroup,SingleGroupQD * SGQD);
    int assertSteadyStateEGoldinP1BC(int iR,int iZ,int iEq,\
        int energyGroup,SingleGroupQD * SGQD);

    // wrapper to assert either a flux or current BC
    void assertSteadyStateWBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertSteadyStateEBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertSteadyStateNBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    void assertSteadyStateSBC(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);

    // functions to enforce coefficients for calculation of steady state 
    // facial currents
    int calcSteadyStateSouthCurrent(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int calcSteadyStateNorthCurrent(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int calcSteadyStateWestCurrent(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);
    int calcSteadyStateEastCurrent(int iR,int iZ,int iEq,int energyGroup,\
        SingleGroupQD * SGQD);

    int steadyStateGreyGroupSources(int iR,int iZ,int iEq,int toEnergyGroup,\
        vector<double> geoParams);

  private:
    // private variables
    YAML::Node * input;
    Mesh * mesh;
    Materials * materials;

    // indices for accessing index and geometry parameters vectors
    const int iCF = 0;
    const int iWF = 1, iEF = 2, iNF = 3, iSF = 4;
    const int iWC = 5, iEC = 6, iNC = 7, iSC = 8;
};

//==============================================================================

#endif
