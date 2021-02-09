#ifndef GREYGROUPSOLVER_H
#define GREYGROUPSOLVER_H 

#include "Mesh.h"
#include "Materials.h"
#include "MultiPhysicsCoupledQD.h"
#include "MultiGroupDNP.h"

class GreyGroupQD;

using namespace std; 

//==============================================================================
//! GreyGroupSolver class that solves RZ neutron transport

class GreyGroupSolver
{
  public:
    // public functions
    GreyGroupSolver(GreyGroupQD * myGGQD,\
        Mesh * myMesh,\
        Materials * myMaterials,\
        YAML::Node * myInput);
    void formLinearSystem();
    void formBackCalcSystem();
    void formSteadyStateLinearSystem();
    void formSteadyStateBackCalcSystem();

    // functions to map grid indices to global index
    vector<int> getIndices(int iR,int iZ);

    // functions to calculate geometry parameters
    double calcVolAvgR(double rDown,double rUp);
    vector<double> calcGeoParams(int iR,int iZ);

    // functions to get Eddington factors
    double getWestErr(int iZ,int iR);
    double getWestErz(int iZ,int iR);
    double getEastErr(int iZ,int iR);
    double getEastErz(int iZ,int iR);
    double getNorthEzz(int iZ,int iR);
    double getNorthErz(int iZ,int iR);
    double getSouthEzz(int iZ,int iR);
    double getSouthErz(int iZ,int iR);

    // functions to enforce governing equations
    void assertZerothMoment(int iR,int iZ,int iEq);
    void applyRadialBoundary(int iR,int iZ,int iEq);
    void applyAxialBoundary(int iR,int iZ,int iEq);

    // functions to enforce governing equations
    void assertSteadyStateZerothMoment(int iR,int iZ,int iEq);
    void applySteadyStateRadialBoundary(int iR,int iZ,int iEq);
    void applySteadyStateAxialBoundary(int iR,int iZ,int iEq);

    // functions to enforce coefficients for facial currents
    void southCurrent(double coeff,int iR,int iZ,int iEq);
    void northCurrent(double coeff,int iR,int iZ,int iEq);
    void westCurrent(double coeff,int iR,int iZ,int iEq);
    void eastCurrent(double coeff,int iR,int iZ,int iEq);

    // functions to enforce coefficients for steady state facial currents
    void steadyStateSouthCurrent(double coeff,int iR,int iZ,int iEq);
    void steadyStateNorthCurrent(double coeff,int iR,int iZ,int iEq);
    void steadyStateWestCurrent(double coeff,int iR,int iZ,int iEq);
    void steadyStateEastCurrent(double coeff,int iR,int iZ,int iEq);

    // functions to enforce coefficients for calculation of facial currents
    void calcSouthCurrent(int iR,int iZ,int iEq);
    void calcNorthCurrent(int iR,int iZ,int iEq);
    void calcWestCurrent(int iR,int iZ,int iEq);
    void calcEastCurrent(int iR,int iZ,int iEq);

    // functions to enforce coefficients for calculation of facial currents
    void calcSteadyStateSouthCurrent(int iR,int iZ,int iEq);
    void calcSteadyStateNorthCurrent(int iR,int iZ,int iEq);
    void calcSteadyStateWestCurrent(int iR,int iZ,int iEq);
    void calcSteadyStateEastCurrent(int iR,int iZ,int iEq);

    // functions to assert flux BCs
    void assertWFluxBC(int iR,int iZ,int iEq);
    void assertEFluxBC(int iR,int iZ,int iEq);
    void assertNFluxBC(int iR,int iZ,int iEq);
    void assertSFluxBC(int iR,int iZ,int iEq);

    // functions to assert current BCs
    void assertWCurrentBC(int iR,int iZ,int iEq);
    void assertECurrentBC(int iR,int iZ,int iEq);
    void assertNCurrentBC(int iR,int iZ,int iEq);
    void assertSCurrentBC(int iR,int iZ,int iEq);

    // functions to assert steady state current BCs
    void assertSteadyStateWCurrentBC(int iR,int iZ,int iEq);
    void assertSteadyStateECurrentBC(int iR,int iZ,int iEq);
    void assertSteadyStateNCurrentBC(int iR,int iZ,int iEq);
    void assertSteadyStateSCurrentBC(int iR,int iZ,int iEq);

    // functions to assert Gol'din BCs
    void assertNGoldinBC(int iR,int iZ,int iEq);
    void assertSGoldinBC(int iR,int iZ,int iEq);
    void assertEGoldinBC(int iR,int iZ,int iEq);

    // functions to assert Gol'din diffusion BCs
    void assertNGoldinP1BC(int iR,int iZ,int iEq);
    void assertSGoldinP1BC(int iR,int iZ,int iEq);
    void assertEGoldinP1BC(int iR,int iZ,int iEq);
    
    // functions to assert steady state Gol'din BCs
    void assertSteadyStateNGoldinBC(int iR,int iZ,int iEq);
    void assertSteadyStateSGoldinBC(int iR,int iZ,int iEq);
    void assertSteadyStateEGoldinBC(int iR,int iZ,int iEq);

    // functions to assert steady state Gol'din BCs
    void assertSteadyStateNGoldinP1BC(int iR,int iZ,int iEq);
    void assertSteadyStateSGoldinP1BC(int iR,int iZ,int iEq);
    void assertSteadyStateEGoldinP1BC(int iR,int iZ,int iEq);

    // wrapper to assert either a flux or current BC
    void assertWBC(int iR,int iZ,int iEq);
    void assertEBC(int iR,int iZ,int iEq);
    void assertNBC(int iR,int iZ,int iEq);
    void assertSBC(int iR,int iZ,int iEq);
    
    // wrapper to assert either a steady state flux or current BC
    void assertSteadyStateWBC(int iR,int iZ,int iEq);
    void assertSteadyStateEBC(int iR,int iZ,int iEq);
    void assertSteadyStateNBC(int iR,int iZ,int iEq);
    void assertSteadyStateSBC(int iR,int iZ,int iEq);

    double calcScatterAndFissionCoeff(int iR,int iZ);
    double calcIntegratingFactor(int iR,int iZ,double rEval,int iLoc);

    // function to solve linear system
    void backCalculateCurrent();

    // check optional inputs
    void checkOptionalParams();

    // function to parse solution vector
    int getFlux();
    int getCurrent();
    int setFlux();
    Eigen::VectorXd getFluxSolutionVector();
    Eigen::VectorXd getCurrentSolutionVector();

    // function to assign pointers 
    void assignPointers(Eigen::SparseMatrix<double,Eigen::RowMajor> * myA,\
        Eigen::VectorXd * myx,\
        Eigen::VectorXd * myxpast,\
        Eigen::VectorXd * myb);

    // public variables
    Eigen::SparseMatrix<double,Eigen::RowMajor> * A;
    Eigen::VectorXd * b;
    Eigen::VectorXd * x;
    Eigen::VectorXd * xPast;
    Eigen::SparseMatrix<double> C;
    Eigen::SparseMatrix<double,Eigen::RowMajor> Atemp;
    Eigen::VectorXd xFlux;
    Eigen::VectorXd currPast;
    Eigen::VectorXd d;
    int energyGroups,nR,nZ,nUnknowns,nCurrentUnknowns;
    bool reflectingBCs = false;
    bool goldinBCs = false;
    bool diffusionBCs = false;

    /* PETSc variables and functions */
    // Variables
    Mat A_p; 
    Vec x_p,xPast_p,b_p;
    Vec xPast_p_seq;
    Mat C_p;
    Vec currPast_p,d_p,xFlux_p;
    Vec currPast_p_seq;
    KSP ksp;
    PC pc;

    // Functions

    void assignMPQDPointer(MultiPhysicsCoupledQD * myMPQD);

    /* STEADY STATE FUNCTIONS */

    void formSteadyStateLinearSystem_p();
    int formSteadyStateBackCalcSystem_p();
    int backCalculateCurrent_p();

    // // functions to enforce governing equations
    int assertSteadyStateZerothMoment_p(int iR,int iZ,int iEq);
    void applySteadyStateRadialBoundary_p(int iR,int iZ,int iEq);
    void applySteadyStateAxialBoundary_p(int iR,int iZ,int iEq);

    // functions to enforce coefficients for steady state facial currents
    int steadyStateSouthCurrent_p(double coeff,int iR,int iZ,int iEq);
    int steadyStateNorthCurrent_p(double coeff,int iR,int iZ,int iEq);
    int steadyStateWestCurrent_p(double coeff,int iR,int iZ,int iEq);
    int steadyStateEastCurrent_p(double coeff,int iR,int iZ,int iEq);
    
    // functions to enforce coefficients for calculation of facial currents
    int calcSteadyStateSouthCurrent_p(int iR,int iZ,int iEq);
    int calcSteadyStateNorthCurrent_p(int iR,int iZ,int iEq);
    int calcSteadyStateWestCurrent_p(int iR,int iZ,int iEq);
    int calcSteadyStateEastCurrent_p(int iR,int iZ,int iEq);

    // wrapper to assert either a steady state flux or current BC
    void assertSteadyStateWBC_p(int iR,int iZ,int iEq);
    void assertSteadyStateEBC_p(int iR,int iZ,int iEq);
    void assertSteadyStateNBC_p(int iR,int iZ,int iEq);
    void assertSteadyStateSBC_p(int iR,int iZ,int iEq);
    
    // functions to assert steady state current BCs
    void assertSteadyStateWCurrentBC_p(int iR,int iZ,int iEq);
    void assertSteadyStateECurrentBC_p(int iR,int iZ,int iEq);
    void assertSteadyStateNCurrentBC_p(int iR,int iZ,int iEq);
    void assertSteadyStateSCurrentBC_p(int iR,int iZ,int iEq);

    // functions to assert steady state Gol'din BCs
    int assertSteadyStateNGoldinBC_p(int iR,int iZ,int iEq);
    int assertSteadyStateSGoldinBC_p(int iR,int iZ,int iEq);
    int assertSteadyStateEGoldinBC_p(int iR,int iZ,int iEq);

    // functions to assert steady state Gol'din BCs
    int assertSteadyStateNGoldinP1BC_p(int iR,int iZ,int iEq);
    int assertSteadyStateSGoldinP1BC_p(int iR,int iZ,int iEq);
    int assertSteadyStateEGoldinP1BC_p(int iR,int iZ,int iEq);

    // functions to assert flux BCs
    int assertWFluxBC_p(int iR,int iZ,int iEq);
    int assertEFluxBC_p(int iR,int iZ,int iEq);
    int assertNFluxBC_p(int iR,int iZ,int iEq);
    int assertSFluxBC_p(int iR,int iZ,int iEq);

    /* TRANSIENT FUNCTIONS */
    
    void formLinearSystem_p();
    int formBackCalcSystem_p();

    // functions to enforce governing equations
    int assertZerothMoment_p(int iR,int iZ,int iEq);
    void applyRadialBoundary_p(int iR,int iZ,int iEq);
    void applyAxialBoundary_p(int iR,int iZ,int iEq);

    // functions to enforce coefficients for facial currents
    int southCurrent_p(double coeff,int iR,int iZ,int iEq);
    int northCurrent_p(double coeff,int iR,int iZ,int iEq);
    int westCurrent_p(double coeff,int iR,int iZ,int iEq);
    int eastCurrent_p(double coeff,int iR,int iZ,int iEq);

    // functions to enforce coefficients for calculation of facial currents
    int calcSouthCurrent_p(int iR,int iZ,int iEq);
    int calcNorthCurrent_p(int iR,int iZ,int iEq);
    int calcWestCurrent_p(int iR,int iZ,int iEq);
    int calcEastCurrent_p(int iR,int iZ,int iEq);
    
    // wrapper to assert either a flux or current BC
    void assertWBC_p(int iR,int iZ,int iEq);
    void assertEBC_p(int iR,int iZ,int iEq);
    void assertNBC_p(int iR,int iZ,int iEq);
    void assertSBC_p(int iR,int iZ,int iEq);

    // functions to assert current BCs
    void assertWCurrentBC_p(int iR,int iZ,int iEq);
    void assertECurrentBC_p(int iR,int iZ,int iEq);
    void assertNCurrentBC_p(int iR,int iZ,int iEq);
    void assertSCurrentBC_p(int iR,int iZ,int iEq);

    // functions to assert Gol'din BCs
    int assertNGoldinBC_p(int iR,int iZ,int iEq);
    int assertSGoldinBC_p(int iR,int iZ,int iEq);
    int assertEGoldinBC_p(int iR,int iZ,int iEq);

    // functions to assert Gol'din diffusion BCs
    int assertNGoldinP1BC_p(int iR,int iZ,int iEq);
    int assertSGoldinP1BC_p(int iR,int iZ,int iEq);
    int assertEGoldinP1BC_p(int iR,int iZ,int iEq);

  private:
    // private variables
    GreyGroupQD * GGQD;
    MultiPhysicsCoupledQD * MPQD;
    YAML::Node * input;
    Mesh * mesh;
    Materials * materials;
    // indices for accessing index and geometry parameters vectors

    int iCF = 0;
    int iWF = 1, iEF = 2, iNF = 3, iSF = 4;
    int iWC = 5, iEC = 6, iNC = 7, iSC = 8;
};

//==============================================================================

#endif
