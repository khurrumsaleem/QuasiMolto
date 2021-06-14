#ifndef MultiGroupDNP_H
#define MultiGroupDNP_H

#include "Mesh.h"
#include "Materials.h"
#include "PETScWrapper.h"

using namespace std;

class MultiPhysicsCoupledQD;
class SingleGroupDNP; 

//==============================================================================
//! Container for all precursor group

class MultiGroupDNP
{
  public:
    int indexOffset = 0;
    int nCoreUnknowns,nRecircUnknowns;
    Eigen::VectorXd beta;
    vector< shared_ptr<SingleGroupDNP> > DNPs; 
    Eigen::VectorXd recircb,recircx;
    Eigen::SparseMatrix<double,Eigen::RowMajor> recircA;
    Eigen::MatrixXd dnpSource;
    MultiGroupDNP(Materials * myMats,\
        Mesh * myMesh,\
        YAML::Node * myInput,\
        MultiPhysicsCoupledQD * myMPQD,\
        int myIndexOffset);
    void readInput();
    void buildCoreLinearSystem();
    void buildSteadyStateCoreLinearSystem();
    void buildRecircLinearSystem();
    void buildSteadyStateRecircLinearSystem();
    void getCoreDNPConc();
    void getCumulativeDNPDecaySource();
    void setCoreDNPConc();
    void printCoreDNPConc();
    void getRecircDNPConc();
    void setRecircDNPConc();
    void printRecircDNPConc();
    void solveRecircLinearSystem();
    MultiPhysicsCoupledQD * mpqd;
    YAML::Node * input;

    /* PETSc */
    Mat recircA_p;
    Vec recircx_p,recircb_p;
    KSP ksp;
    PC pc;
   
    // Dual purpose
    int solveRecircLinearSystem_p();

    // Steady state
    void buildSteadyStateCoreLinearSystem_p();
    int buildSteadyStateRecircLinearSystem_p();
    
    // Transient
    void buildCoreLinearSystem_p();
    int buildRecircLinearSystem_p();

    // Pseudo transient
    void buildPseudoTransientCoreLinearSystem_p();

  private:
    Materials * mats;
    Mesh * mesh; 

};

//==============================================================================

#endif
