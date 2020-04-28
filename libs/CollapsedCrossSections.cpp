/// File: CollapsedCrossSections.cpp
// Purpose: contain group collapsed nuclear data
// Date: October 28, 2019

#include "CollapsedCrossSections.h"

using namespace std; 
using namespace arma;

///==============================================================================
/// CollapsedCrossSections class object constructor
///
/// @param [in] myMats Materials object for this simulation
CollapsedCrossSections::CollapsedCrossSections(int nZ,int nR)
{
  
  // Initialize matrices holding one-group cross sections
  sigT.setZero(nZ,nR);  
  sigS.setZero(nZ,nR);  
  sigF.setZero(nZ,nR);  
  chiP.setZero(nZ,nR);  
  chiD.setZero(nZ,nR);  
  neutV.setZero(nZ,nR);  
  chiD.setZero(nZ,nR);  

};
//==============================================================================

