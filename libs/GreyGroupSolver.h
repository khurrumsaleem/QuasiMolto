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

    // functions to map grid indices to global index
    vector<int> getIndices(int iR,int iZ);

    // functions to calculate geometry parameters
    double calcVolAvgR(double rDown,double rUp);
    vector<double> calcGeoParams(int iR,int iZ);

    double calcScatterAndFissionCoeff(int iR,int iZ);
    double calcIntegratingFactor(int iR,int iZ,double rEval,int iLoc);

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

    // functions to get Eddington factors
    double getWestErr(int iZ,int iR);
    double getWestErz(int iZ,int iR);
    double getEastErr(int iZ,int iR);
    double getEastErz(int iZ,int iR);
    double getNorthEzz(int iZ,int iR);
    double getNorthErz(int iZ,int iR);
    double getSouthEzz(int iZ,int iR);
    double getSouthErz(int iZ,int iR);

    /* STEADY STATE FUNCTIONS */

    void formSteadyStateLinearSystem();
    int formSteadyStateBackCalcSystem();
    int backCalculateCurrent();

    // // functions to enforce governing equations
    int assertSteadyStateZerothMoment(int iR,int iZ,int iEq);
    void applySteadyStateRadialBoundary(int iR,int iZ,int iEq);
    void applySteadyStateAxialBoundary(int iR,int iZ,int iEq);

    // functions to enforce coefficients for steady state facial currents
    int steadyStateSouthCurrent(double coeff,int iR,int iZ,int iEq);
    int steadyStateNorthCurrent(double coeff,int iR,int iZ,int iEq);
    int steadyStateWestCurrent(double coeff,int iR,int iZ,int iEq);
    int steadyStateEastCurrent(double coeff,int iR,int iZ,int iEq);
    
    // functions to enforce coefficients for calculation of facial currents
    int calcSteadyStateSouthCurrent(int iR,int iZ,int iEq);
    int calcSteadyStateNorthCurrent(int iR,int iZ,int iEq);
    int calcSteadyStateWestCurrent(int iR,int iZ,int iEq);
    int calcSteadyStateEastCurrent(int iR,int iZ,int iEq);

    // wrapper to assert either a steady state flux or current BC
    void assertSteadyStateWBC(int iR,int iZ,int iEq);
    void assertSteadyStateEBC(int iR,int iZ,int iEq);
    void assertSteadyStateNBC(int iR,int iZ,int iEq);
    void assertSteadyStateSBC(int iR,int iZ,int iEq);
    
    // functions to assert steady state current BCs
    void assertSteadyStateWCurrentBC(int iR,int iZ,int iEq);
    void assertSteadyStateECurrentBC(int iR,int iZ,int iEq);
    void assertSteadyStateNCurrentBC(int iR,int iZ,int iEq);
    void assertSteadyStateSCurrentBC(int iR,int iZ,int iEq);

    // functions to assert steady state Gol'din BCs
    int assertSteadyStateNGoldinBC(int iR,int iZ,int iEq);
    int assertSteadyStateSGoldinBC(int iR,int iZ,int iEq);
    int assertSteadyStateEGoldinBC(int iR,int iZ,int iEq);

    // functions to assert steady state Gol'din BCs
    int assertSteadyStateNGoldinP1BC(int iR,int iZ,int iEq);
    int assertSteadyStateSGoldinP1BC(int iR,int iZ,int iEq);
    int assertSteadyStateEGoldinP1BC(int iR,int iZ,int iEq);

    // functions to assert flux BCs
    int assertWFluxBC(int iR,int iZ,int iEq);
    int assertEFluxBC(int iR,int iZ,int iEq);
    int assertNFluxBC(int iR,int iZ,int iEq);
    int assertSFluxBC(int iR,int iZ,int iEq);

    /* TRANSIENT FUNCTIONS */
    
    void formLinearSystem();
    int formBackCalcSystem();

    // functions to enforce governing equations
    int assertZerothMoment(int iR,int iZ,int iEq);
    void applyRadialBoundary(int iR,int iZ,int iEq);
    void applyAxialBoundary(int iR,int iZ,int iEq);

    // functions to enforce coefficients for facial currents
    int southCurrent(double coeff,int iR,int iZ,int iEq);
    int northCurrent(double coeff,int iR,int iZ,int iEq);
    int westCurrent(double coeff,int iR,int iZ,int iEq);
    int eastCurrent(double coeff,int iR,int iZ,int iEq);

    // functions to enforce coefficients for calculation of facial currents
    int calcSouthCurrent(int iR,int iZ,int iEq);
    int calcNorthCurrent(int iR,int iZ,int iEq);
    int calcWestCurrent(int iR,int iZ,int iEq);
    int calcEastCurrent(int iR,int iZ,int iEq);
    
    // wrapper to assert either a flux or current BC
    void assertWBC(int iR,int iZ,int iEq);
    void assertEBC(int iR,int iZ,int iEq);
    void assertNBC(int iR,int iZ,int iEq);
    void assertSBC(int iR,int iZ,int iEq);

    // functions to assert current BCs
    void assertWCurrentBC(int iR,int iZ,int iEq);
    void assertECurrentBC(int iR,int iZ,int iEq);
    void assertNCurrentBC(int iR,int iZ,int iEq);
    void assertSCurrentBC(int iR,int iZ,int iEq);

    // functions to assert Gol'din BCs
    int assertNGoldinBC(int iR,int iZ,int iEq);
    int assertSGoldinBC(int iR,int iZ,int iEq);
    int assertEGoldinBC(int iR,int iZ,int iEq);

    // functions to assert Gol'din diffusion BCs
    int assertNGoldinP1BC(int iR,int iZ,int iEq);
    int assertSGoldinP1BC(int iR,int iZ,int iEq);
    int assertEGoldinP1BC(int iR,int iZ,int iEq);

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
