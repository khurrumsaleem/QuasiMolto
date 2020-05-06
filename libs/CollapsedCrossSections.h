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
    Eigen::MatrixXd sigT,sigS,sigF,chiP,chiD,neutV,nu,qdFluxCoeff; 
    vector<Eigen::MatrixXd> groupDNPFluxCoeff; 
  
    // Functions 
    double dnpFluxCoeff(int iZ,int iR,int dnpID);
};

//==============================================================================

#endif
