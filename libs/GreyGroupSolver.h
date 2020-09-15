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
    
    // functions to assert steady state Gol'din BCs
    void assertSteadyStateNGoldinBC(int iR,int iZ,int iEq);
    void assertSteadyStateSGoldinBC(int iR,int iZ,int iEq);
    void assertSteadyStateEGoldinBC(int iR,int iZ,int iEq);

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
    void getFlux();
    void getCurrent();
    void setFlux();
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

  private:
    // private variables
    GreyGroupQD * GGQD;
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
