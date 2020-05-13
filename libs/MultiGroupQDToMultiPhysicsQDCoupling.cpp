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

  // Store past values
  mats->oneGroupXS->zNeutVPast = mats->oneGroupXS->zNeutV;
  mats->oneGroupXS->rNeutVPast = mats->oneGroupXS->rNeutV;

  // Calculate flux and current weighted data
  calculateFluxWeightedData();
  calculateFluxWeightedBCData();
  calculateAxialCurrentWeightedData();
  calculateRadialCurrentWeightedData();

  // Calculate zeta factors
  calculateAxialZetaFactors(); 
  calculateRadialZetaFactors(); 
};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux weighting 
///
void MGQDToMPQDCoupling::calculateFluxWeightedData()
{

  // Temporary accumulator variables
  double fluxAccum,totalSigS;
  double flux,beta,nu;
  double mySigT,mySigS,mySigF,myNeutV,myErr,myEzz,myErz;

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
        mats->oneGroupXS->sigT(iZ,iR) += mySigT*flux;

        // Flux weighted sigS
        mats->oneGroupXS->sigS(iZ,iR) += mySigS*flux;

        // Flux weighted neutron velocity
        mats->oneGroupXS->neutV(iZ,iR) += myNeutV*flux;

        // Flux weighted Ezz
        mats->oneGroupXS->Ezz(iZ,iR) += myEzz*flux;

        // Flux weighted Err
        mats->oneGroupXS->Err(iZ,iR) += myErr*flux;

        // Flux weighted Erz
        mats->oneGroupXS->Erz(iZ,iR) += myErz*flux;

        // Flux weighted quasidiffusion coefficient
        mats->oneGroupXS->qdFluxCoeff(iZ,iR) += (1-beta)*nu*mySigF*flux;

        // Loop over DNP groups
        for (int iDNPGroup = 0; iDNPGroup < mpqd->mgdnp->DNPs.size(); 
            iDNPGroup++) 
        {
          // Flux weighted DNP coefficient 
          beta = mpqd->mgdnp->DNPs[iDNPGroup]->beta(iEnergyGroup);
          mats->oneGroupXS->groupDNPFluxCoeff[iDNPGroup](iZ,iR)\
            += beta*nu*mySigF*flux; 
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

};
//==============================================================================

//==============================================================================
/// Collapse boundary condition data 
///
void MGQDToMPQDCoupling::calculateFluxWeightedBCData()
{

  // Temporary accumulator variables
  double nFluxAccum,nFlux,nInwardFlux,nRatio;
  double sFluxAccum,sFlux,sInwardFlux,sRatio;
  double eFluxAccum,eFlux,eInwardFlux,eRatio;
  double nInwardCurrent,sInwardCurrent,eInwardCurrent;

  // Reset ratios
  mpqd->ggqd->nOutwardCurrToFluxRatioBC.setZero();
  mpqd->ggqd->sOutwardCurrToFluxRatioBC.setZero();
  mpqd->ggqd->eOutwardCurrToFluxRatioBC.setZero();

  mpqd->ggqd->nOutwardCurrToFluxRatioInwardWeightedBC.setZero();
  mpqd->ggqd->sOutwardCurrToFluxRatioInwardWeightedBC.setZero();
  mpqd->ggqd->eOutwardCurrToFluxRatioInwardWeightedBC.setZero();

  mpqd->ggqd->nInwardFluxBC.setZero();
  mpqd->ggqd->sInwardFluxBC.setZero();
  mpqd->ggqd->eInwardFluxBC.setZero();

  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR; iR++)
  {

    // Reset accumulators
    nFluxAccum = 0.0;
    sFluxAccum = 0.0;

    // Loop over neutron energy groups
    for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
    {

      // Get flux in these cells and energy group
      nFlux = mgqd->SGQDs[iEnergyGroup]->sFluxZ(0,iR); 
      sFlux = mgqd->SGQDs[iEnergyGroup]->sFluxZ(mesh->nZ,iR); 

      nInwardFlux = mgqd->SGQDs[iEnergyGroup]->nInwardFluxBC(iR); 
      sInwardFlux = mgqd->SGQDs[iEnergyGroup]->sInwardFluxBC(iR); 

      nInwardCurrent = mgqd->SGQDs[iEnergyGroup]->nInwardCurrentBC(iR); 
      sInwardCurrent = mgqd->SGQDs[iEnergyGroup]->sInwardCurrentBC(iR); 

      // Get outward current to flux ratio for these cells
      nRatio = mgqd->SGQDs[iEnergyGroup]->nOutwardCurrToFluxRatioBC(iR); 
      sRatio = mgqd->SGQDs[iEnergyGroup]->sOutwardCurrToFluxRatioBC(iR); 

      // Flux weighted quasidiffusion coefficient
      mpqd->ggqd->nOutwardCurrToFluxRatioBC(iR) += nRatio*nFlux; 
      mpqd->ggqd->sOutwardCurrToFluxRatioBC(iR) += sRatio*sFlux;

      mpqd->ggqd->nOutwardCurrToFluxRatioInwardWeightedBC(iR)\
        += nRatio*nInwardFlux; 
      mpqd->ggqd->sOutwardCurrToFluxRatioInwardWeightedBC(iR)\
        += sRatio*sInwardFlux;

      // Accumulate fluxes and currents 
      nFluxAccum = nFluxAccum + nFlux;
      sFluxAccum = sFluxAccum + sFlux;

      mpqd->ggqd->nInwardFluxBC(iR) += nInwardFlux;
      mpqd->ggqd->sInwardFluxBC(iR) += sInwardFlux;

      mpqd->ggqd->nInwardCurrentBC(iR) += nInwardCurrent;
      mpqd->ggqd->sInwardCurrentBC(iR) += sInwardCurrent;

    }

    // Divide data by accumulated values
    mpqd->ggqd->nOutwardCurrToFluxRatioBC(iR)\
      = mpqd->ggqd->nOutwardCurrToFluxRatioBC(iR)/nFluxAccum; 
    mpqd->ggqd->sOutwardCurrToFluxRatioBC(iR)\
      = mpqd->ggqd->sOutwardCurrToFluxRatioBC(iR)/sFluxAccum; 

    mpqd->ggqd->nOutwardCurrToFluxRatioInwardWeightedBC(iR)\
      = mpqd->ggqd->nOutwardCurrToFluxRatioInwardWeightedBC(iR)\
      /mpqd->ggqd->nInwardFluxBC(iR); 
    mpqd->ggqd->sOutwardCurrToFluxRatioInwardWeightedBC(iR)\
      = mpqd->ggqd->sOutwardCurrToFluxRatioInwardWeightedBC(iR)\
      /mpqd->ggqd->sInwardFluxBC(iR); 

  }

  // Loop over spatial mesh
  for (int iZ = 0; iZ < mesh->nZ; iZ++)
  {

    // Reset accumulators
    eFluxAccum = 0.0;

    // Loop over neutron energy groups
    for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
    {

      // Get flux in these cells and energy group
      eFlux = mgqd->SGQDs[iEnergyGroup]->sFluxR(iZ,mesh->nR); 

      eInwardFlux = mgqd->SGQDs[iEnergyGroup]->eInwardFluxBC(iZ); 

      eInwardCurrent = mgqd->SGQDs[iEnergyGroup]->eInwardCurrentBC(iZ); 

      // Get outward current to flux ratio for these cells
      eRatio = mgqd->SGQDs[iEnergyGroup]->eOutwardCurrToFluxRatioBC(iZ); 

      // Flux weighted quasidiffusion coefficient
      mpqd->ggqd->eOutwardCurrToFluxRatioBC(iZ) += eRatio*nFlux;

      mpqd->ggqd->eOutwardCurrToFluxRatioInwardWeightedBC(iZ)\
        += eRatio*nInwardFlux;

      // Accumulate fluxes and currents 
      eFluxAccum = eFluxAccum + eFlux;

      mpqd->ggqd->eInwardFluxBC(iZ) += eInwardFlux;

      mpqd->ggqd->eInwardCurrentBC(iZ) += eInwardCurrent;
    }

    // Divide data by accumulated values
    mpqd->ggqd->eOutwardCurrToFluxRatioBC(iZ)\
      = mpqd->ggqd->eOutwardCurrToFluxRatioBC(iZ)/nFluxAccum; 

    mpqd->ggqd->eOutwardCurrToFluxRatioInwardWeightedBC(iZ)\
      = mpqd->ggqd->eOutwardCurrToFluxRatioInwardWeightedBC(iZ)\
      /mpqd->ggqd->eInwardFluxBC(iZ); 
  }

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with axial current weighting 
///
void MGQDToMPQDCoupling::calculateAxialCurrentWeightedData()
{

  // Temporary accumulator variables
  double zCurrentAccum;
  double zCurrent;
  double mySigT,myNeutV;

  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ+1; iZ++)
    {

      // Reset accumulator
      zCurrentAccum = 0.0;

      // Loop over neutron energy groups
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        zCurrent = abs(mgqd->SGQDs[iEnergyGroup]->currentZ(iZ,iR));

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
        mats->oneGroupXS->zSigTR(iZ,iR) += mySigT*zCurrent;

        // Axial current weighted sigTR
        mats->oneGroupXS->zNeutV(iZ,iR) += zCurrent/myNeutV;

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
};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with radial current weighting 
///
void MGQDToMPQDCoupling::calculateRadialCurrentWeightedData()
{

  // Temporary accumulator variables
  double rCurrentAccum;
  double rCurrent;
  double mySigT,myNeutV;

  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR+1; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ; iZ++)
    {

      // Reset accumulators
      rCurrentAccum = 0.0;

      // Loop over neutron energy groups
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        rCurrent = abs(mgqd->SGQDs[iEnergyGroup]->currentR(iZ,iR)); 

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
        mats->oneGroupXS->rSigTR(iZ,iR) += mySigT*rCurrent;

        // Radial current weighted neutron velocity
        mats->oneGroupXS->rNeutV(iZ,iR) += rCurrent/myNeutV;

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

//==============================================================================
/// Calculate radial zeta factors 
///
void MGQDToMPQDCoupling::calculateRadialZetaFactors()
{

  // Temporary accumulator variables
  double fluxAccum;
  double flux,rCurrent,rCurrentPast;
  double mySigT,mySigTR,myNeutV,myRNeutV,myRNeutVPast;
  double timeDerivative,sigTDiff,pastSum,presentSum; 


  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR+1; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ; iZ++)
    {

      // Reset accumulators
      fluxAccum = 0.0;
      pastSum = 0.0;
      presentSum = 0.0;
      sigTDiff = 0.0;

      // Now that the current weight neutron velocities are calculated, we 
      // calculate we can calculate the zeta factors
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        rCurrent = mgqd->SGQDs[iEnergyGroup]->currentR(iZ,iR); 
        rCurrentPast = mgqd->SGQDs[iEnergyGroup]->currentRPrev(iZ,iR); 
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
        mySigTR = mats->oneGroupXS->rSigTR(iZ,iR);
        myNeutV = mats->neutV(iEnergyGroup);
        myRNeutV = mats->oneGroupXS->rNeutV(iZ,iR);
        myRNeutVPast = mats->oneGroupXS->rNeutVPast(iZ,iR);

        // Radial current weighted zeta1
        pastSum = pastSum + (1.0/myNeutV-1.0/myRNeutVPast)*rCurrentPast;
        presentSum = presentSum + (1.0/myNeutV-1.0/myRNeutV)*rCurrent;

        // Radial current weighted zeta2
        sigTDiff = sigTDiff + (mySigT-mySigTR)*rCurrent;       

        // Accumulate absolute radial current
        fluxAccum = fluxAccum + flux;

      }

      // Approximate derivative
      timeDerivative = (presentSum-pastSum)/mesh->dt; 

      // Divide by accumulated flux
      mats->oneGroupXS->rZeta1(iZ,iR) = timeDerivative/fluxAccum; 
      mats->oneGroupXS->rZeta2(iZ,iR) = sigTDiff/fluxAccum; 
      mats->oneGroupXS->rZeta(iZ,iR) = mats->oneGroupXS->rZeta1(iZ,iR)\
                                       + mats->oneGroupXS->rZeta2(iZ,iR); 
    }
  }
};
//==============================================================================

//==============================================================================
/// Calculate axial zeta factors 
///
void MGQDToMPQDCoupling::calculateAxialZetaFactors()
{

  // Temporary accumulator variables
  double fluxAccum;
  double flux,zCurrent,zCurrentPast;
  double mySigT,mySigTR,myNeutV,myZNeutV,myZNeutVPast;
  double timeDerivative,sigTDiff,pastSum,presentSum; 


  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ+1; iZ++)
    {

      // Reset accumulators
      fluxAccum = 0.0;
      pastSum = 0.0;
      presentSum = 0.0;
      sigTDiff = 0.0;

      // Now that the current weight neutron velocities are calculated, we 
      // calculate we can calculate the zeta factors
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        zCurrent = mgqd->SGQDs[iEnergyGroup]->currentZ(iZ,iR); 
        zCurrentPast = mgqd->SGQDs[iEnergyGroup]->currentZPrev(iZ,iR); 
        flux = mgqd->SGQDs[iEnergyGroup]->sFluxZ(iZ,iR);

        // Get nuclear data at these locations and energy groups
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
        mySigTR = mats->oneGroupXS->zSigTR(iZ,iR);
        myNeutV = mats->neutV(iEnergyGroup);
        myZNeutV = mats->oneGroupXS->zNeutV(iZ,iR);
        myZNeutVPast = mats->oneGroupXS->zNeutVPast(iZ,iR);

        // Radial current weighted zeta1
        pastSum = pastSum + (1.0/myNeutV-1.0/myZNeutVPast)*zCurrentPast;
        presentSum = presentSum + (1.0/myNeutV-1.0/myZNeutV)*zCurrent;

        // Radial current weighted zeta2
        sigTDiff = sigTDiff + (mySigT-mySigTR)*zCurrent;       

        // Accumulate absolute radial current
        fluxAccum = fluxAccum + flux;

      }

      // Approximate derivative
      timeDerivative = (presentSum-pastSum)/mesh->dt; 

      // Divide by accumulated flux
      mats->oneGroupXS->zZeta1(iZ,iR) = timeDerivative/fluxAccum; 
      mats->oneGroupXS->zZeta2(iZ,iR) = sigTDiff/fluxAccum; 
      mats->oneGroupXS->zZeta(iZ,iR) = mats->oneGroupXS->zZeta1(iZ,iR)\
                                       + mats->oneGroupXS->zZeta2(iZ,iR); 
    }
  }
};
//==============================================================================

