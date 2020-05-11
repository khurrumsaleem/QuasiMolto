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
  double flux,rCurrent,zCurrent,beta,nu;
  double mySigT,mySigS,mySigF,myNeutV,myErr,myEzz,myErz;

  // Reset collapsed nuclear data
  mats->oneGroupXS->resetData();

  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ; iZ++)
    {

      // Reset accumulators
      fluxAccum = 0.0;

      // Loop over neutron energy groups
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        flux = mgqd->SGQDs[iEnergyGroup]->sFlux(iZ,iR); 
      
        // Get beta and nu for this energy group 
        beta = mpqd->mgdnp->beta(iEnergyGroup); 
        nu = mats->nu(iZ,iR,iEnergyGroup); 

        // Get nuclear data at these locations and energy groups
        mySigT = mats->sigT(iZ,iR,iEnergyGroup);
        mySigS= 0; 
        for (int iScat = 0; iScat < mats->nGroups; iScat++)
        {
          mySigS = totalSigS + mats->sigS(iZ,iR,0,iScat); 
        }
        mySigF = mats->sigF(iZ,iR,iEnergyGroup);
        myNeutV = mats->neutV(iEnergyGroup);
        myErr = mgqd->SGQDs[iEnergyGroup]->Err(iZ,iR);
        myEzz = mgqd->SGQDs[iEnergyGroup]->Ezz(iZ,iR);
        myErz = mgqd->SGQDs[iEnergyGroup]->Erz(iZ,iR);

        // Flux weighted sigT
        mats->oneGroupXS->sigT(iZ,iR) = mySigT*flux;

        // Flux weighted sigS
        mats->oneGroupXS->sigS(iZ,iR) = mySigS*flux;

        // Flux weighted neutron velocity
        mats->oneGroupXS->neutV(iZ,iR) = myNeutV*flux;

        // Flux weighted Ezz
        mats->oneGroupXS->Ezz(iZ,iR) = myEzz*flux;

        // Flux weighted Err
        mats->oneGroupXS->Err(iZ,iR) = myErr*flux;

        // Flux weighted Erz
        mats->oneGroupXS->Erz(iZ,iR) = myErz*flux;

        // Flux weighted quasidiffusion coefficient
        mats->oneGroupXS->qdFluxCoeff(iZ,iR) = (1-beta)*nu*mySigF*flux;

        // Loop over DNP groups
        for (int iDNPGroup = 0; iDNPGroup < mpqd->mgdnp->DNPs.size(); 
            iDNPGroup++) 
        {
          // Flux weighted DNP coefficient 
          beta = mpqd->mgdnp->DNPs[iDNPGroup]->beta(iEnergyGroup);
          mats->oneGroupXS->groupDNPFluxCoeff[iDNPGroup](iZ,iR)\
            = beta*nu*mySigF*flux; 
        }

        // Accumulate flux
        fluxAccum = fluxAccum + flux;
        
      }

      // Divide data by accumulated values
      mats->oneGroupXS->sigT(iZ,iR) = mats->oneGroupXS->sigT(iZ,iR)/fluxAccum;
      mats->oneGroupXS->sigS(iZ,iR) = mats->oneGroupXS->sigS(iZ,iR)/fluxAccum;
      mats->oneGroupXS->neutV(iZ,iR) = fluxAccum/mats->oneGroupXS->neutV(iZ,iR);
      mats->oneGroupXS->Ezz(iZ,iR) = mats->oneGroupXS->Ezz(iZ,iR)/fluxAccum;
      mats->oneGroupXS->Err(iZ,iR) = mats->oneGroupXS->Err(iZ,iR)/fluxAccum;
      mats->oneGroupXS->Erz(iZ,iR) = mats->oneGroupXS->Erz(iZ,iR)/fluxAccum;
      mats->oneGroupXS->qdFluxCoeff(iZ,iR) = mats->oneGroupXS->\
                                             qdFluxCoeff(iZ,iR)/fluxAccum;

      // Loop over DNP groups
      for (int iDNPGroup = 0; iDNPGroup < mpqd->mgdnp->DNPs.size(); 
          iDNPGroup++) 
      {
        mats->oneGroupXS->groupDNPFluxCoeff[iDNPGroup](iZ,iR)\
          = mats->oneGroupXS->groupDNPFluxCoeff[iDNPGroup](iZ,iR)\
          /fluxAccum; 
      }

    }
  }

  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ+1; iZ++)
    {

      // Reset accumulator
      zCurrentAccum = 0.0;
      fluxAccum = 0.0;

      // Loop over neutron energy groups
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        zCurrent = abs(mgqd->SGQDs[iEnergyGroup]->currentZ(iZ,iR));
        flux = mgqd->SGQDs[iEnergyGroup]->sFluxZ(iZ,iR);

        // Get nuclear data at these locations and energy groups
        // Assuming a matching material exist on the other side of boundaries 
        if (iZ == 0)
          mySigT = mats->sigT(iZ,iR,iEnergyGroup);
        else if (iZ == mesh->nZ)
          mySigT = mats->sigT(iZ-1,iR,iEnergyGroup);
        else
        {
          mySigT = (mats->sigT(iZ-1,iR,iEnergyGroup)*mesh->dzsCorner(iZ-1)\
              + mats->sigT(iZ,iR,iEnergyGroup)*mesh->dzsCorner(iZ))\
                   /(mesh->dzsCorner(iZ-1)+mesh->dzsCorner(iZ));
        }     
        myNeutV = mats->neutV(iEnergyGroup);

        // Axial current weighted neutron velocity
        mats->oneGroupXS->zSigTR(iZ,iR) = mySigT*zCurrent;

        // Axial current weighted sigTR
        mats->oneGroupXS->zNeutV(iEnergyGroup) = myNeutV*zCurrent;

        // Axial current weighted zeta1
        // Axial current weighted zeta2

        // Accumulate absolute axial current
        zCurrentAccum = zCurrentAccum + zCurrent;
        
      }
      // Divide data by accumulated values
      mats->oneGroupXS->zSigTR(iZ,iR) = mats->oneGroupXS->zSigTR(iZ,iR)\
                                        /zCurrentAccum;
      mats->oneGroupXS->zNeutV(iZ,iR) = zCurrentAccum\
                                        /mats->oneGroupXS->zNeutV(iZ,iR);
    }
  }

  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR+1; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ; iZ++)
    {

      // Reset accumulators
      rCurrentAccum = 0.0;
      fluxAccum = 0.0;

      // Loop over neutron energy groups
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        rCurrent = abs(mgqd->SGQDs[iEnergyGroup]->currentR(iZ,iR)); 
        flux = mgqd->SGQDs[iEnergyGroup]->sFluxR(iZ,iR);
      
        // Get nuclear data at these locations and energy groups
        if (iR == 0)
          mySigT = mats->sigT(iZ,iR,iEnergyGroup);
        else if (iR == mesh->nR)
          mySigT = mats->sigT(iZ,iR-1,iEnergyGroup);
        else
        {
          mySigT = (mats->sigT(iZ,iR-1,iEnergyGroup)*mesh->drsCorner(iR-1)\
              + mats->sigT(iZ,iR,iEnergyGroup)*mesh->drsCorner(iR))\
                   /(mesh->drsCorner(iR-1)+mesh->drsCorner(iR));
        }     
        myNeutV = mats->neutV(iEnergyGroup);

        // Radial current weighted sigTR
        mats->oneGroupXS->rSigTR(iZ,iR) = mySigT*rCurrent;

        // Radial current weighted neutron velocity
        mats->oneGroupXS->rNeutV(iZ,iR) = myNeutV*rCurrent;

        // Radial current weighted zeta1
        // Radial current weighted zeta2

        // Accumulate absolute radial current
        rCurrentAccum = rCurrentAccum + rCurrent;
        
      }
      // Divide data by accumulated values
      mats->oneGroupXS->rSigTR(iZ,iR) = mats->oneGroupXS->rSigTR(iZ,iR)\
                                        /rCurrentAccum;
      mats->oneGroupXS->rNeutV(iZ,iR) = rCurrentAccum\
                                        /mats->oneGroupXS->rNeutV(iZ,iR);
      
    }
  }

};
//==============================================================================
