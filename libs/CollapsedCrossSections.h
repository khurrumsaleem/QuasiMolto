#ifndef COLLAPSEDCROSSSECTIONS_H
#define COLLAPSEDCROSSSECTIONS_H

#include "Mesh.h"

using namespace std; 
using namespace arma;

//==============================================================================
//! Class the holds collapsed one-group nuclear data

class CollapsedCrossSections
{
  public:
    
    // Constructor
    CollapsedCrossSections(int nZ,int nR);
    
    // Variables
    Eigen::MatrixXd sigT,sigS,sigF,rSigTR,zSigTR,neutV,rNeutV,\
      rNeutVPast,zNeutV,zNeutVPast,qdFluxCoeff,Ezz,Err,Erz,rZeta1,rZeta2,\
      rZeta,zZeta1,zZeta2,zZeta; 
    vector<Eigen::MatrixXd> groupDNPFluxCoeff; 
  
    // Functions 
    double dnpFluxCoeff(int iZ,int iR,int dnpID);
    void print();
    void resetData();
};

//==============================================================================

#endif
