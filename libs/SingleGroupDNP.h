#ifndef SingleGroupDNP_H
#define SingleGroupDNP_H

#include "Mesh.h"
#include "Materials.h"

using namespace std;

class MultiGroupDNP; // forward declaration
class GreyGroupQD; // forward declaration

//==============================================================================
//! Contains information and builds linear system for a single precursor group

class SingleGroupDNP
{
  public:
    Eigen::SparseMatrix<double,Eigen::RowMajor> Atemp;
    Eigen::MatrixXd dnpConc,recircConc,flux,recircFlux,dirac,recircDirac;
    Eigen::MatrixXd inletConc,recircInletConc;
    Eigen::VectorXd inletVelocity,recircInletVelocity,outletConc,recircOutletConc;
    Eigen::VectorXd beta; 
    double lambda;
    int coreInletIndex,coreOutletIndex,recircInletIndex,recircOutletIndex;
    string fluxLimiter = "superbee";
    int coreIndexOffset = 0;
    int recircIndexOffset = 0;
    int dnpID = 0;
    SingleGroupDNP(Materials * myMats,\
        Mesh * myMesh,\
        MultiGroupDNP * myMGDNPS,\
        Eigen::VectorXd myBeta,\
        double myLambda,\
        int myCoreIndexOffset,\
        int myRecircIndexOffset,\
        int myDNPid);
    int getIndex(int iZ,int iR,int indexOffset);
    void assignBoundaryIndices();
    void updateBoundaryConditions();
    void calcCoreDNPFluxes();
    void calcRecircDNPFluxes();
    void buildCoreLinearSystem();
    void buildSteadyStateCoreLinearSystem();
    void buildRecircLinearSystem();
    void buildSteadyStateRecircLinearSystem();
    void buildLinearSystem(Eigen::SparseMatrix<double,Eigen::RowMajor> * myA,\
        Eigen::VectorXd * myb,\
        Eigen::MatrixXd myDNPConc,\
        Eigen::MatrixXd myDNPFlux,\
        arma::rowvec dzs,\
        int myIndexOffset,
        bool fluxSource = true);
    void buildSteadyStateLinearSystem(\
        Eigen::SparseMatrix<double,Eigen::RowMajor> * myA,\
        Eigen::VectorXd * myb,\
        Eigen::MatrixXd myDNPConc,\
        Eigen::MatrixXd myDNPFlux,\
        Eigen::MatrixXd myInletDNP,\
        arma::rowvec dzs,\
        int myIndexOffset,
        bool fluxSource = true);

    Eigen::MatrixXd getInitialConc(double initConc);
    Eigen::MatrixXd calcDiracs(Eigen::MatrixXd dnpConc,\
        Eigen::MatrixXd inletConc,\
        Eigen::VectorXd outletConc);
    Eigen::MatrixXd calcFluxes(Eigen::MatrixXd myDNPConc,\
        Eigen::MatrixXd myFlowVelocity,\
        Eigen::MatrixXd dirac,\
        Eigen::MatrixXd inletConc,\
        Eigen::VectorXd inletVelocity,\
        arma::rowvec dzs);
    Eigen::MatrixXd calcImplicitFluxes(Eigen::MatrixXd myDNPConc,\
        Eigen::MatrixXd myFlowVelocity,\
        Eigen::MatrixXd inletConc,\
        Eigen::VectorXd inletVelocity,\
        arma::rowvec dzs);
    int getCoreConc();
    void setCoreConc();
    int getRecircConc();
    void setRecircConc();
    double calcPhi(double theta,string fluxLimiter); 
    double calcTheta(double DNPupwindInterface,double DNPinterface);

    /* PETSc functions */   

    // Steady state

    int buildSteadyStateCoreLinearSystem_p();
    int buildSteadyStateRecircLinearSystem_p();
    int buildSteadyStateLinearSystem_p(\
        Mat * myA_p,\
        Vec * myb_p,\
        Eigen::MatrixXd myDNPConc,\
        Eigen::MatrixXd myDNPFlux,\
        Eigen::MatrixXd myInletDNP,\
        arma::rowvec dzs,\
        int myIndexOffset,
        bool fluxSource = true);

  private: 
    Materials * mats;
    Mesh * mesh;
    MultiGroupDNP * mgdnp;
};

//==============================================================================

#endif
