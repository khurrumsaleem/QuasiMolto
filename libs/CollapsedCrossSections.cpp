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
CollapsedCrossSections::CollapsedCrossSections(int nZ,int nR,int nEnergyGroups)
{
  
  // Initialize matrices holding one-group cross sections
  sigT.setZero(nZ,nR);  
  sigS.setZero(nZ,nR);  
  sigF.setZero(nZ,nR);  
  rSigTR.setZero(nZ,nR+1);  
  zSigTR.setZero(nZ+1,nR);  
  neutV.setZero(nZ,nR);  
  rNeutV.setZero(nZ,nR+1);  
  rNeutVPast.setZero(nZ,nR+1);  
  zNeutV.setZero(nZ+1,nR);  
  zNeutVPast.setZero(nZ+1,nR);  
  Ezz.setConstant(nZ,nR,1.0/3.0);  
  Err.setConstant(nZ,nR,1.0/3.0);  
  Erz.setZero(nZ,nR);  
  rZeta1.setZero(nZ,nR+1);  
  rZeta2.setZero(nZ,nR+1);  
  zZeta1.setZero(nZ+1,nR);  
  zZeta2.setZero(nZ+1,nR);  
  rZeta.setZero(nZ,nR+1);  
  zZeta.setZero(nZ+1,nR);  
  qdFluxCoeff.setZero(nZ,nR);  

  // Set size of matrices contained in groupSigS 
  groupSigS.resize(nEnergyGroups);
  for (int iGroup = 0; iGroup < groupSigS.size(); iGroup++)
    groupSigS[iGroup].setZero(nZ,nR);

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
/// Return scattering cross section into neutron energy group eID at location

/// @param [in] nZ number of axial cells
/// @param [in] nR number of radial cells
/// @param [in] eID energy group index 
double CollapsedCrossSections::groupScatterXS(int iZ,int iR,int eIdx)
{

  return groupSigS[eIdx](iZ,iR);

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
  zZeta.setZero();  
  rZeta.setZero();  
  qdFluxCoeff.setZero();  

  // Reset dnpFluxCoeff
  for (int iDNPGroup = 0; iDNPGroup < groupDNPFluxCoeff.size(); iDNPGroup++)
    groupDNPFluxCoeff[iDNPGroup].setZero();
 
  // Reset group scattering cross sections 
  for (int iGroup = 0; iGroup < groupSigS.size(); iGroup++)
    groupSigS[iGroup].setZero();

};
//==============================================================================

///==============================================================================
/// Print collapsed nuclear data 
///
void CollapsedCrossSections::print()
{

  // Set Eigen::MatrixXd types to zeros
  cout << "sigT: " << endl;
  cout << sigT  << endl;
  cout << endl;

  cout << "sigS: " << endl;
  cout << sigS  << endl;
  cout << endl;

  cout << "sigF: " << endl;
  cout << sigF  << endl;
  cout << endl;

  cout << "rSigTR: " << endl;
  cout << rSigTR  << endl;
  cout << endl;

  cout << "zSigTR: " << endl;
  cout << zSigTR  << endl;
  cout << endl;

  cout << "neutV: " << endl;
  cout << neutV  << endl;
  cout << endl;

  cout << "rNeutV: " << endl;
  cout << rNeutV  << endl;
  cout << endl;

  cout << "zNeutV: " << endl;
  cout << zNeutV  << endl;
  cout << endl;

  cout << "rNeutVPast: " << endl;
  cout << rNeutVPast  << endl;
  cout << endl;

  cout << "zNeutVPast: " << endl;
  cout << zNeutVPast  << endl;
  cout << endl;

  cout << "Ezz: " << endl;
  cout << Ezz  << endl;
  cout << endl;

  cout << "Err: " << endl;
  cout << Err  << endl;
  cout << endl;

  cout << "Erz: " << endl;
  cout << Erz  << endl;
  cout << endl;

  cout << "rZeta1: " << endl;
  cout << rZeta1  << endl;
  cout << endl;

  cout << "rZeta2: " << endl;
  cout << rZeta2  << endl;
  cout << endl;

  cout << "rZeta: " << endl;
  cout << rZeta  << endl;
  cout << endl;

  cout << "zZeta1: " << endl;
  cout << zZeta1  << endl;
  cout << endl;

  cout << "zZeta2: " << endl;
  cout << zZeta2  << endl;
  cout << endl;

  cout << "zZeta: " << endl;
  cout << zZeta  << endl;
  cout << endl;

  cout << "qdFluxCoeff: " << endl;
  cout << qdFluxCoeff  << endl;
  cout << endl;

  // Print dnpFluxCoeff
  cout << "groupDNPFluxCoeff: " << endl;
  for (int iDNPGroup = 0; iDNPGroup < groupDNPFluxCoeff.size(); iDNPGroup++)
  {
    cout << groupDNPFluxCoeff[iDNPGroup] << endl;
    cout << endl;
  }

};
//==============================================================================
