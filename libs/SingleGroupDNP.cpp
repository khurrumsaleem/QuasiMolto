// File: SingleGroupDNP.cpp     
// Purpose: Contains information and builds linear system for a precursor group
// Date: April 9, 2020

#include "SingleGroupDNP.h"
#include "MultiGroupDNP.h"
#include "MultiPhysicsCoupledQD.h"

using namespace std;

//==============================================================================
/// SingleGroupDNP class object constructor
///
/// @param [in] blankType blank for this material
SingleGroupDNP::SingleGroupDNP(Materials * myMats,\
  Mesh * myMesh,\
  MultiGroupDNP * myMGDNP,\
  double myBeta,\
  double myLambda)
{

  // Assign inputs to their member variables
  mats = myMats;
  mesh = myMesh;
  mgdnp = myMGDNP;
  beta = myBeta;
  lambda = myLambda;  
};
//==============================================================================

//==============================================================================
/// Calculate diracs to model advection of precursors
///
void SingleGroupDNP::calcDiracs()
{
};
//==============================================================================

//==============================================================================
/// Calculate fluxes to model advection of precursors
///
void SingleGroupDNP::calcFluxes()
{
};
//==============================================================================


//==============================================================================
/// Map 2D coordinates to index of delayed neutron precursor concentration in 
/// the 1D solution vector
///
/// @param [in] iZ axial location
/// @param [in] iR radial location
/// @param [out] index the index for dnp concentration in the 1D solution vector
int SingleGroupDNP::getIndex(int iZ, int iR)
{

  int index,nR = mesh->drsCorner.size();

  index = indexOffset + iR + nR*iZ;

  return index;
  
};
//==============================================================================

//==============================================================================
/// Map 2D coordinates to index of delayed neutron precursor concentration in 
/// the 1D solution vector
///
/// @param [in] iZ axial location
/// @param [in] iR radial location
/// @param [out] index the index for temperature in the 1D solution vector
void SingleGroupDNP::getConc()
{

  for (int iR = 0; iR < mesh-> drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh-> dzsCorner.size(); iZ++)
    { 
    
      dnpConc(iZ,iR) = mgdnp->mpqd->x(getIndex(iZ,iR));   
    
    }
  }
  
};
//==============================================================================
