#ifndef COLLAPSEDCROSSSECTIONS_H
#define COLLAPSEDCROSSSECTIONS_H

#include "Mesh.h"
#include "WriteData.h"

using namespace std; 
using namespace arma;

//==============================================================================
//! Class the holds collapsed one-group nuclear data

class CollapsedCrossSections
{
  public:
    
    // Constructor
    CollapsedCrossSections(Mesh * myMesh,int nEnergyGroups);
    
    // Variables
    // sigS is the flux weighted one group scattering cross section
    // groupSigS is the flux weighted one group to group scattering cross 
    // section
    Eigen::MatrixXd sigT,sigS,sigF,rSigTR,zSigTR,neutV,neutVPast,rNeutV,\
      rNeutVPast,zNeutV,zNeutVPast,qdFluxCoeff,Ezz,Err,Erz,rZeta1,rZeta2,\
      rZeta,zZeta1,zZeta2,zZeta; 
    vector<Eigen::MatrixXd> groupDNPFluxCoeff; 
    vector<Eigen::MatrixXd> groupSigS; 
    vector<Eigen::MatrixXd> groupUpscatterCoeff; 
    string outputDir = "1GXS/";
    Mesh * mesh;
  
    // Functions 
    double dnpFluxCoeff(int iZ,int iR,int dnpID);
    double groupScatterXS(int iZ,int iR,int dnpID);
    double upscatterCoeff(int iZ,int iR,int dnpID);
    void print();
    void writeVars();
    void resetData();
};

//==============================================================================

#endif
