#ifndef GREYGROUPSOLVERBASE_H
#define GREYGROUPSOLVERBASE_H 

#include "Mesh.h"
#include "Materials.h"
#include "MultiPhysicsCoupledQD.h"
#include "MultiGroupDNP.h"

class GreyGroupQD;

using namespace std; 

//==============================================================================
//! GreyGroupSolver class that solves RZ neutron transport

class GreyGroupSolverBase
{
  public:

    // Default constructor
    GreyGroupSolverBase(GreyGroupQD * myGGQD,\
        Mesh * myMesh,\
        Materials * myMaterials,\
        YAML::Node * myInput);

    GreyGroupSolverBase(const GreyGroupSolverBase &solver);

    // Public variables
    int energyGroups,nR,nZ,nUnknowns,nCurrentUnknowns;
    bool reflectingBCs = false;
    bool goldinBCs = false;
    bool diffusionBCs = false;

    /* PETSc variables and functions */
    Mat * A; 
    Vec * x,*xPast,*b;
    Vec * xPastSeq;
    Mat C;
    Vec currPast,d,xFlux;
    Vec currPastSeq;
    KSP ksp;
    PC pc;
    Eigen::SparseMatrix<double,Eigen::RowMajor> Atemp;

    // Public functions
    void formLinearSystem();
    int formBackCalcSystem();
    int backCalculateCurrent();
    void checkOptionalParams();
    void assignMPQDPointer(MultiPhysicsCoupledQD * myMPQD);
    void assignLinearSystemPointers(Mat * myA,\
        Vec * myx,\
        Vec * myxPast,\
        Vec * myxPastSeq,\
        Vec * myb);
    
    // funcations to apply interface conditions
    void applyRadialBoundary(int iR,int iZ,int iEq);
    void applyAxialBoundary(int iR,int iZ,int iEq);

    // functions to assert flux BCs
    int assertWFluxBC(int iR,int iZ,int iEq);
    int assertEFluxBC(int iR,int iZ,int iEq);
    int assertNFluxBC(int iR,int iZ,int iEq);
    int assertSFluxBC(int iR,int iZ,int iEq);

    // functions to assert current BCs
    void assertWCurrentBC(int iR,int iZ,int iEq);
    void assertECurrentBC(int iR,int iZ,int iEq);
    void assertNCurrentBC(int iR,int iZ,int iEq);
    void assertSCurrentBC(int iR,int iZ,int iEq);

    // functions to assert Gol'din BCs
    int assertNGoldinBC(int iR,int iZ,int iEq);
    int assertSGoldinBC(int iR,int iZ,int iEq);
    int assertEGoldinBC(int iR,int iZ,int iEq);

    // functions to assert Gol'din P1 BCs
    int assertNGoldinP1BC(int iR,int iZ,int iEq);
    int assertSGoldinP1BC(int iR,int iZ,int iEq);
    int assertEGoldinP1BC(int iR,int iZ,int iEq);
    
    // wrapper to assert either a flux or current BC
    void assertWBC(int iR,int iZ,int iEq);
    void assertEBC(int iR,int iZ,int iEq);
    void assertNBC(int iR,int iZ,int iEq);
    void assertSBC(int iR,int iZ,int iEq);

    int getFlux();
    int getCurrent();
    int setFlux();
    vector<int> getIndices(int iR,int iZ);

    double getWestErr(int iZ,int iR);
    double getWestErz(int iZ,int iR);
    double getEastErr(int iZ,int iR);
    double getEastErz(int iZ,int iR);
    double getNorthEzz(int iZ,int iR);
    double getNorthErz(int iZ,int iR);
    double getSouthEzz(int iZ,int iR);
    double getSouthErz(int iZ,int iR);
    double calcVolAvgR(double rDown,double rUp);
    double calcScatterAndFissionCoeff(int iR,int iZ);
    double calcIntegratingFactor(int iR,int iZ,double rEval,int iLoc);
    vector<double> calcGeoParams(int iR,int iZ);

    Eigen::VectorXd getFluxSolutionVector();
    Eigen::VectorXd getCurrentSolutionVector();

  protected:

    // Private variables
    GreyGroupQD * GGQD;
    MultiPhysicsCoupledQD * MPQD;
    YAML::Node * input;
    Mesh * mesh;
    Materials * materials;

    // indices for accessing index and geometry parameters vectors
    int iCF = 0;
    int iWF = 1, iEF = 2, iNF = 3, iSF = 4;
    int iWC = 5, iEC = 6, iNC = 7, iSC = 8;
    
    // functions to enforce governing equations
    virtual int assertZerothMoment(int iR,int iZ,int iEq) = 0;

    // functions to enforce coefficients for facial currents
    virtual int southCurrent(double coeff,int iR,int iZ,int iEq) = 0;
    virtual int northCurrent(double coeff,int iR,int iZ,int iEq) = 0;
    virtual int westCurrent(double coeff,int iR,int iZ,int iEq) = 0;
    virtual int eastCurrent(double coeff,int iR,int iZ,int iEq) = 0;

    // functions to enforce coefficients for calculation of facial currents
    virtual int calcSouthCurrent(int iR,int iZ,int iEq) = 0;
    virtual int calcNorthCurrent(int iR,int iZ,int iEq) = 0;
    virtual int calcWestCurrent(int iR,int iZ,int iEq) = 0;
    virtual int calcEastCurrent(int iR,int iZ,int iEq) = 0;

};

//==============================================================================

#endif
