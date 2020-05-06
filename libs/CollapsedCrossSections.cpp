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
  chiP.setZero(nZ,nR);  
  chiD.setZero(nZ,nR);  
  neutV.setZero(nZ,nR);  
  chiD.setZero(nZ,nR);  
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
