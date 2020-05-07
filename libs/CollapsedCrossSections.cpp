/// File: CollapsedCrossSections.cpp
// Purpose: contain group collapsed nuclear data
// Date: October 28, 2019

#include "CollapsedCrossSections.h"

using namespace std; 
using namespace arma;

///==============================================================================
/// CollapsedCrossSections class object constructor
///
/// @param [in] nZ number of axial cells
/// @param [in] nR number of radial cells
CollapsedCrossSections::CollapsedCrossSections(int nZ,int nR)
{
  
  // Initialize matrices holding one-group cross sections
  sigT.setZero(nZ,nR);  
  sigS.setZero(nZ,nR);  
  sigF.setZero(nZ,nR);  
  rSigTR.setZero(nZ,nR+1);  
  zSigTR.setZero(nZ+1,nR);  
  chiP.setZero(nZ,nR);  
  chiD.setZero(nZ,nR);  
  neutV.setZero(nZ,nR);  
  rNeutV.setZero(nZ,nR+1);  
  zNeutV.setZero(nZ+1,nR);  
  Ezz.setZero(nZ,nR);  
  Err.setZero(nZ,nR);  
  Erz.setZero(nZ,nR);  
  rZeta1.setZero(nZ,nR+1);  
  rZeta2.setZero(nZ,nR+1);  
  zZeta1.setZero(nZ+1,nR);  
  zZeta2.setZero(nZ+1,nR);  
  nu.setZero(nZ,nR);  
  qdFluxCoeff.setZero(nZ,nR);  

};
//==============================================================================

///==============================================================================
/// Return DNP flux coefficient at a location and precursor group
///
/// @param [in] nZ number of axial cells
/// @param [in] nR number of radial cells
/// @param [in] dnpID number of precursor group
double CollapsedCrossSections::dnpFluxCoeff(int iZ,int iR,int dnpID)
{

  return groupDNPFluxCoeff[dnpID](iZ,iR);

};
//==============================================================================
