#ifndef SINGLEGROUPQD_H 
#define SINGLEGROUPQD_H

#include "Mesh.h"
#include "Materials.h"
#include "MultiGroupQD.h"
#include "QuasidiffusionSolver.h"

using namespace std;
using namespace arma;

//==============================================================================
//! SingleGroupQD class that holds quasidiffusion information

class SingleGroupQD
{
  public:
    int energyGroup;
    Eigen::MatrixXd sFlux;
    Eigen::MatrixXd sFluxR;
    Eigen::MatrixXd sFluxRPrev;
    Eigen::MatrixXd sFluxZ;
    Eigen::MatrixXd sFluxZPrev;
    Eigen::MatrixXd currentR;
    Eigen::MatrixXd currentRPrev;
    Eigen::MatrixXd currentZ;
    Eigen::MatrixXd currentZPrev;
    Eigen::MatrixXd sFluxPrev;
    Eigen::MatrixXd q;
    Eigen::MatrixXd fissionSource;
    Eigen::MatrixXd scatterSource;

    // Eddington factors
    Eigen::MatrixXd Err,ErrPrev;
    Eigen::MatrixXd Ezz,EzzPrev;
    Eigen::MatrixXd Erz,ErzPrev;

    // flux boundary conditions  
    Eigen::VectorXd wFluxBC;
    Eigen::VectorXd eFluxBC;
    Eigen::VectorXd nFluxBC;
    Eigen::VectorXd sFluxBC;

    // current boundary conditions
    Eigen::VectorXd wCurrentRBC;
    Eigen::VectorXd eCurrentRBC;
    Eigen::VectorXd nCurrentZBC;
    Eigen::VectorXd sCurrentZBC;

    // vectors for robin boundary conditions
    Eigen::VectorXd eInwardCurrentBC;
    Eigen::VectorXd nInwardCurrentBC;
    Eigen::VectorXd sInwardCurrentBC;

    Eigen::VectorXd eInwardFluxBC;
    Eigen::VectorXd nInwardFluxBC;
    Eigen::VectorXd sInwardFluxBC;

    Eigen::VectorXd eOutwardCurrToFluxRatioBC;
    Eigen::VectorXd nOutwardCurrToFluxRatioBC;
    Eigen::VectorXd sOutwardCurrToFluxRatioBC;

    Eigen::VectorXd eAbsCurrentBC;
    Eigen::VectorXd nAbsCurrentBC;
    Eigen::VectorXd sAbsCurrentBC;

    // public functions
    SingleGroupQD(int myEnergyGroup,\
        MultiGroupQD * myMGQD,\
        Materials * myMaterials,\
        Mesh * myMesh,\
        YAML::Node * myInput);
    void formContributionToLinearSystem();
    void formContributionToBackCalcSystem();
    void getFlux();
    Eigen::VectorXd getFluxSolutionVector();
    Eigen::VectorXd getCurrentSolutionVector();
    void checkOptionalParams();
    void writeFlux();

  private:
    MultiGroupQD * MGQD;
    Materials * mats;
    YAML::Node * input;
    Mesh * mesh;

};

//==============================================================================

#endif
