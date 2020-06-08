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
CollapsedCrossSections::CollapsedCrossSections(Mesh * myMesh, int nEnergyGroups)
{

  // Assign pointer 
  mesh = myMesh;

  // Declare initialization variables 
  int nZ = mesh->nZ;
  int nR = mesh->nR;

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
  groupUpscatterCoeff.resize(nEnergyGroups);
  for (int iGroup = 0; iGroup < groupSigS.size(); iGroup++)
  {
    groupSigS[iGroup].setZero(nZ,nR);
    groupUpscatterCoeff[iGroup].setZero(nZ,nR);
  }

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
/// @param [in] eIdx energy group index 
double CollapsedCrossSections::groupScatterXS(int iZ,int iR,int eIdx)
{

  return groupSigS[eIdx](iZ,iR);

};
//==============================================================================

///==============================================================================
/// Return scattering cross section into neutron energy group eID at location

/// @param [in] nZ number of axial cells
/// @param [in] nR number of radial cells
/// @param [in] eIdx energy group index 
double CollapsedCrossSections::upscatterCoeff(int iZ,int iR,int eIdx)
{

  return groupUpscatterCoeff[eIdx](iZ,iR);

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
  {
    groupSigS[iGroup].setZero();
    groupUpscatterCoeff[iGroup].setZero();
  }

};
//==============================================================================

///==============================================================================
/// Write variables to output directory 
///
void CollapsedCrossSections::writeVars()
{
  
  string varName;

  // Write all simple 1GXS variables 
  mesh->output->write(outputDir,"sigT",sigT);

  mesh->output->write(outputDir,"sigS",sigS);

  mesh->output->write(outputDir,"sigF",sigF);

  mesh->output->write(outputDir,"rSigTR",rSigTR);

  mesh->output->write(outputDir,"zSigTR",zSigTR);

  mesh->output->write(outputDir,"neutV",neutV);

  mesh->output->write(outputDir,"rNeutV",rNeutV);

  mesh->output->write(outputDir,"zNeutV",zNeutV);

  mesh->output->write(outputDir,"rZeta1",rZeta1);

  mesh->output->write(outputDir,"rZeta2",rZeta2);

  mesh->output->write(outputDir,"rZeta",rZeta);

  mesh->output->write(outputDir,"zZeta1",rZeta1);

  mesh->output->write(outputDir,"zZeta2",rZeta2);

  mesh->output->write(outputDir,"zZeta",rZeta);

  mesh->output->write(outputDir,"qdFluxCoeff",qdFluxCoeff);

  // Write dnpFluxCoeff
  for (int iDNPGroup = 0; iDNPGroup < groupDNPFluxCoeff.size(); iDNPGroup++)
  {
    varName = "DNPFluxCoeff_Group_" + to_string(iDNPGroup);
    mesh->output->write(outputDir,varName,groupDNPFluxCoeff[iDNPGroup]);
  }

  // Write groupSigS 
  for (int iScatter= 0; iScatter < groupSigS.size(); iScatter++)
  {
    varName = "SigS_Group_" + to_string(iScatter);
    mesh->output->write(outputDir,varName,groupSigS[iScatter]);
  }

  // Write groupUpscatterCoeff 
  for (int iScatter= 0; iScatter < groupUpscatterCoeff.size(); iScatter++)
  {
    varName = "UpscatterCoeff_Group__" + to_string(iScatter);
    mesh->output->write(outputDir,varName,groupUpscatterCoeff[iScatter]);
  }

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

  // Print groupSigS 
  cout << "groupSigS: " << endl;
  for (int iScatter= 0; iScatter < groupSigS.size(); iScatter++)
  {
    cout << groupSigS[iScatter] << endl;
    cout << endl;
  }

  // Print groupUpscatterCoeff 
  cout << "groupUpscatterCoeff: " << endl;
  for (int iScatter= 0; iScatter < groupUpscatterCoeff.size(); iScatter++)
  {
    cout << groupUpscatterCoeff[iScatter] << endl;
    cout << endl;
  }

};
//==============================================================================
