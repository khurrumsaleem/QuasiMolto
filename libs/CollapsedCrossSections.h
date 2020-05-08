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
    Eigen::MatrixXd sigT,sigS,sigF,rSigTR,zSigTR,chiP,chiD,neutV,rNeutV,zNeutV,\
      nu,qdFluxCoeff,Ezz,Err,Erz,rZeta1,rZeta2,zZeta1,zZeta2; 
    vector<Eigen::MatrixXd> groupDNPFluxCoeff; 
  
    // Functions 
    double dnpFluxCoeff(int iZ,int iR,int dnpID);
    void resetData();
};

//==============================================================================

#endif
