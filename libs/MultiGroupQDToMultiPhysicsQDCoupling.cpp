// File: MultiGroupQDToMultiPhysicsQDCoupling.cpp     
// Purpose: Calculate coupling factors between multi-group quasidiffusion and
//   multi-physics quasidiffusion.
// Date: April 9, 2020

#include "MultiGroupQDToMultiPhysicsQDCoupling.h"

using namespace std;

//==============================================================================
/// MultiGroupQDToMultiPhysicsQDCoupling class object constructor
///
/// @param [in] myMesh mesh object 
/// @param [in] myMats materials object 
/// @param [in] myInput input object 
/// @param [in] myMPQD MultiPhysicsCoupledQD object 
/// @param [in] myMGPQD multigroup quasidiffusion object 
MGQDToMPQDCoupling::MGQDToMPQDCoupling(Mesh * myMesh,\
        Materials * myMats,\
        YAML::Node * myInput,\
        MultiPhysicsCoupledQD * myMPQD,\
        MultiGroupQD * myMGQD)
{

  // Assign inputs to their member variables
  mesh = myMesh; 
  mats = myMats;
  input = myInput;
  mpqd = myMPQD;
  mgqd = myMGQD;

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
void MGQDToMPQDCoupling::collapseNuclearData()
{

  // Temporary accumulator variables
  double fluxAccum,rCurrentAccum,zCurrentAccum,totalSigS;
  double flux,rCurrent,zCurrent;

  // Reset collapsed nuclear data
  mats->oneGroupXS->resetData();

  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ; iZ++)
    {

      // Reset accumulators
      fluxAccum = 0.0;
      rCurrentAccum = 0.0;
      zCurrentAccum = 0.0;

      // Loop over neutron energy groups
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        flux = mgqd->SGQDs[iEnergyGroup]->sFlux(iZ,iR); 
        rCurrent = mgqd->SGQDs[iEnergyGroup]->currentR(iZ,iR); 
        zCurrent = mgqd->SGQDs[iEnergyGroup]->currentZ(iZ,iR); 

        // Flux weighted sigT
        mats->oneGroupXS->sigT(iZ,iR) = mats->sigT(iZ,iR,iEnergyGroup)*flux;

        // Flux weighted sigS
        totalSigS = 0; 
        for (int iScat = 0; iScat < mats->nGroups; iScat++)
        {
          totalSigS = totalSigS + mats->sigS(iZ,iR,0,iScat); 
        }
        mats->oneGroupXS->sigS(iZ,iR) = totalSigS*flux;

        // Flux weighted neutron velocity
        mats->oneGroupXS->neutV(iEnergyGroup) = mats->neutV(iEnergyGroup)*flux;

        // Flux weighted Ezz
        mats->oneGroupXS->Ezz(iZ,iR) = mgqd->SGQDs[iEnergyGroup]->Ezz(iZ,iR)\
                                       *flux;

        // Flux weighted Err
        mats->oneGroupXS->Err(iZ,iR) = mgqd->SGQDs[iEnergyGroup]->Err(iZ,iR)\
                                       *flux;

        // Flux weighted Erz
        mats->oneGroupXS->Erz(iZ,iR) = mgqd->SGQDs[iEnergyGroup]->Erz(iZ,iR)\
                                       *flux;

        // Flux weighted quasidiffusion coefficient
        // Flux weighted DNP coefficient 

        // Axial current weighted neutron velocity
        // Axial current weighted sigTR

        // Radial current weighted neutron velocity
        // Radial current weighted sigTR

        // Radial current weighted zeta1
        // Radial current weighted zeta2

        // Axial current weighted zeta1
        // Axial current weighted zeta2

        // Loop over DNP groups
        for (int iDNPGroup = 0; iDNPGroup < mpqd->mgdnp->DNPs.size(); 
            iDNPGroup++) 
        {
          // Flux weighted DNP coefficient 
        }

        // Accumulate flux
        // Accumulate absolute axial current
        // Accumulate absolute radial current
        
      }
      // Divide data by accumulated values
    }
  }
};
//==============================================================================
