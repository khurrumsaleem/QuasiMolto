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
    void getCoreDNPConc();
    void getCumulativeDNPDecaySource();
    void setCoreDNPConc();
    void printCoreDNPConc();
    void getRecircDNPConc();
    void setRecircDNPConc();
    void printRecircDNPConc();
    MultiPhysicsCoupledQD * mpqd;
    YAML::Node * input;

    /* PETSc */
    Mat recircA_p;
    Vec recircx_p,recircb_p;
    KSP ksp;
    PC pc;
   
    // Dual purpose
    int solveRecircLinearSystem();

    // Steady state
    void buildSteadyStateCoreLinearSystem();
    int buildSteadyStateRecircLinearSystem();
    
    // Transient
    void buildCoreLinearSystem();
    int buildRecircLinearSystem();

    // Pseudo transient
    void buildPseudoTransientCoreLinearSystem();

  private:
    Materials * mats;
    Mesh * mesh; 

};

//==============================================================================

#endif
