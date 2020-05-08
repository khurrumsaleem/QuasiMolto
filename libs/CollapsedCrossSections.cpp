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

///==============================================================================
/// Reset collapsed nuclear data 
///
void CollapsedCrossSections::resetData()
{

  // Set Eigen::MatrixXd types to zeros
  sigT.setZero();  
  sigS.setZero();  
  sigF.setZero();  
  rSigTR.setZero();  
  zSigTR.setZero();  
  chiP.setZero();  
  chiD.setZero();  
  neutV.setZero();  
  rNeutV.setZero();  
  zNeutV.setZero();  
  Ezz.setZero();  
  Err.setZero();  
  Erz.setZero();  
  rZeta1.setZero();  
  rZeta2.setZero();  
  zZeta1.setZero();  
  zZeta2.setZero();  
  nu.setZero();  
  qdFluxCoeff.setZero();  

  // Reset dnpFluxCoeff
  for (int iDNPGroup = 0; iDNPGroup < qdFluxCoeff.size(); iDNPGroup++)
  {
    groupDNPFluxCoeff[iDNPGroup].setZero();
  }

};
//==============================================================================
