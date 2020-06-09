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

  // Assign MGQD pointer in grey group solver
  mgqd->assignMultiPhysicsCoupledQDPointer(mpqd);

  // Assign MGQD pointer in grey group solver
  mpqd->ggqd->assignMGQDPointer(mgqd);

  // Initialize collapsed nuclear data based on initial values stored on MGQD
  // objects 
  initCollapsedNuclearData();

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
bool MGQDToMPQDCoupling::solveOneStep()
{

  Eigen::VectorXd xCurrentIter,xPrevIter,ones,residualVec;
  double residual;

  for (int iStep = 0; iStep < 100; iStep++)
  {
    mgqd->buildLinearSystem();
    mgqd->solveLinearSystem();
    mgqd->buildBackCalcSystem();
    mgqd->backCalculateCurrent();
    mgqd->getFluxes();

    collapseNuclearData();
    xPrevIter = mpqd->x;
    mpqd->buildLinearSystem();
    mpqd->solveLinearSystem();
    xCurrentIter = mpqd->x;
    //cout << "xCurrentIter:" << endl;
    //cout << xCurrentIter << endl;
    residual = calcResidual(xPrevIter,xCurrentIter);
    cout << "          "; 
    cout << "MGQD->MPQD Residual: " << residual <<endl; 
    if (residual < mpqd->epsMPQD)
    {
      return true; 
    }
    mpqd->ggqd->GGSolver->getFlux();
    mpqd->mgdnp->getCumulativeDNPDecaySource();
  }

  return false;

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
void MGQDToMPQDCoupling::solveTransient()
{

  bool converged; 

  for (int iTime = 0; iTime < mesh->dts.size(); iTime++)
  {
    converged = solveOneStep();  
    if (not converged)
    {
      cout << "Multigroup A: " << endl;
      cout << mgqd->QDSolve->A << endl;
      cout << endl;
      cout << "Multigroup b: " << endl;
      cout << mgqd->QDSolve->b << endl;
      cout << endl;
      cout << "Grey group A: " << endl;
      cout << mpqd->A << endl;
      cout << endl;
      cout << "Grey group b: " << endl;
      cout << mpqd->b << endl;
      cout << endl;
      mgqd->updateVarsAfterConvergence(); 
      mpqd->updateVarsAfterConvergence();
      break;
    }
    mgqd->updateVarsAfterConvergence(); 
    mpqd->updateVarsAfterConvergence();
    mpqd->writeVars(); 
    mgqd->writeVars(); 
    mesh->advanceOneTimeStep();
  }

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
void MGQDToMPQDCoupling::initCollapsedNuclearData()
{

  collapseNuclearData();
  mats->oneGroupXS->zNeutVPast = mats->oneGroupXS->zNeutV;
  mats->oneGroupXS->rNeutVPast = mats->oneGroupXS->rNeutV;

};
//==============================================================================


//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
void MGQDToMPQDCoupling::collapseNuclearData()
{

  // Store past values
  mats->oneGroupXS->zNeutVPast = mats->oneGroupXS->zNeutV;
  mats->oneGroupXS->rNeutVPast = mats->oneGroupXS->rNeutV;

  // Reset collapsed nuclear data
  mats->oneGroupXS->resetData();

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
  double fluxAccum;
  double flux,beta,nu;
  double mySigT,mySigS,mySigF,myNeutV,myErr,myEzz,myErz,myGroupSigS;

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
          mySigS += mats->sigS(iZ,iR,iEnergyGroup,iScat); 
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

        // Flux weighted groupSigS
        for (int iScat = 0; iScat < mats->nGroups; iScat++)
        {
          myGroupSigS = mats->sigS(iZ,iR,iEnergyGroup,iScat); 
          mats->oneGroupXS->groupSigS[iScat](iZ,iR) += myGroupSigS*flux;
          if (iScat < iEnergyGroup) 
            mats->oneGroupXS->groupUpscatterCoeff[iScat](iZ,iR)\
              += myGroupSigS*flux;
        }

        // Flux weighted sigF
        mats->oneGroupXS->sigF(iZ,iR) += mySigF*flux;

        // Flux weighted neutron velocity
        mats->oneGroupXS->neutV(iZ,iR) += flux/myNeutV;

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
      for (int iScat = 0; iScat < mats->nGroups; iScat++)
      {
        mats->oneGroupXS->groupSigS[iScat](iZ,iR)\
          = mats->oneGroupXS->groupSigS[iScat](iZ,iR)/fluxAccum;
        mats->oneGroupXS->groupUpscatterCoeff[iScat](iZ,iR)\
          = mats->oneGroupXS->groupUpscatterCoeff[iScat](iZ,iR)/fluxAccum;
      }
      mats->oneGroupXS->sigF(iZ,iR) = mats->oneGroupXS->sigF(iZ,iR)/fluxAccum;
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

  // Set eddington factors in grey group object 
  mpqd->ggqd->Ezz = mats->oneGroupXS->Ezz;
  mpqd->ggqd->Erz = mats->oneGroupXS->Erz;
  mpqd->ggqd->Err = mats->oneGroupXS->Err;

  // Calculate interface Eddington factors
  calculateFluxWeightedInterfaceEddingtons();

};
//==============================================================================

//==============================================================================
/// Collapse interface Eddingtons with flux weighting 
///
void MGQDToMPQDCoupling::calculateFluxWeightedInterfaceEddingtons()
{

  // Temporary accumulator variables
  double fluxAccum;
  double flux;
  double myErr,myEzz,myErz;

  // Reset Eddington factors
  mpqd->ggqd->EzzRadial.setZero();
  mpqd->ggqd->ErrRadial.setZero();
  mpqd->ggqd->ErzRadial.setZero();
  mpqd->ggqd->EzzAxial.setZero();
  mpqd->ggqd->ErrAxial.setZero();
  mpqd->ggqd->ErzAxial.setZero();

  // Loop over radial interface mesh
  for (int iR = 0; iR < mesh->rCornerEdge.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ; iZ++)
    {

      // Reset accumulators
      fluxAccum = 0.0;

      // Loop over neutron energy groups
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        flux = mgqd->SGQDs[iEnergyGroup]->sFluxR(iZ,iR); 

        myErr = mgqd->SGQDs[iEnergyGroup]->ErrRadial(iZ,iR);
        myEzz = mgqd->SGQDs[iEnergyGroup]->EzzRadial(iZ,iR);
        myErz = mgqd->SGQDs[iEnergyGroup]->ErzRadial(iZ,iR);

        // Flux weighted Ezz
        mpqd->ggqd->EzzRadial(iZ,iR) += myEzz*flux;

        // Flux weighted Err
        mpqd->ggqd->ErrRadial(iZ,iR) += myErr*flux;

        // Flux weighted Erz
        mpqd->ggqd->ErzRadial(iZ,iR) += myErz*flux;

        // Accumulate flux
        fluxAccum = fluxAccum + flux;

      }

      mpqd->ggqd->EzzRadial(iZ,iR) = mpqd->ggqd->EzzRadial(iZ,iR)/fluxAccum;
      mpqd->ggqd->ErrRadial(iZ,iR) = mpqd->ggqd->ErrRadial(iZ,iR)/fluxAccum;
      mpqd->ggqd->ErzRadial(iZ,iR) = mpqd->ggqd->ErzRadial(iZ,iR)/fluxAccum;

    }
  }

  // Loop over radial interface mesh
  for (int iR = 0; iR < mesh->nR; iR++)
  {
    for (int iZ = 0; iZ < mesh->zCornerEdge.size(); iZ++)
    {

      // Reset accumulators
      fluxAccum = 0.0;

      // Loop over neutron energy groups
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        flux = mgqd->SGQDs[iEnergyGroup]->sFluxZ(iZ,iR); 

        myErr = mgqd->SGQDs[iEnergyGroup]->ErrAxial(iZ,iR);
        myEzz = mgqd->SGQDs[iEnergyGroup]->EzzAxial(iZ,iR);
        myErz = mgqd->SGQDs[iEnergyGroup]->ErzAxial(iZ,iR);

        // Flux weighted Ezz
        mpqd->ggqd->EzzAxial(iZ,iR) += myEzz*flux;

        // Flux weighted Err
        mpqd->ggqd->ErrAxial(iZ,iR) += myErr*flux;

        // Flux weighted Erz
        mpqd->ggqd->ErzAxial(iZ,iR) += myErz*flux;

        // Accumulate flux
        fluxAccum = fluxAccum + flux;

      }

      mpqd->ggqd->EzzAxial(iZ,iR) = mpqd->ggqd->EzzAxial(iZ,iR)/fluxAccum;
      mpqd->ggqd->ErrAxial(iZ,iR) = mpqd->ggqd->ErrAxial(iZ,iR)/fluxAccum;
      mpqd->ggqd->ErzAxial(iZ,iR) = mpqd->ggqd->ErzAxial(iZ,iR)/fluxAccum;

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
  double nFluxAccum,nFlux,nInwardFlux,nInwardFluxAccum,nRatio;
  double sFluxAccum,sFlux,sInwardFlux,sInwardFluxAccum,sRatio;
  double eFluxAccum,eFlux,eInwardFlux,eInwardFluxAccum,eRatio;
  double nInwardCurrent,sInwardCurrent,eInwardCurrent;
  double nBias,sBias,eBias;
  double nInwardBias,sInwardBias,eInwardBias;

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

  mpqd->ggqd->nInwardCurrentBC.setZero();
  mpqd->ggqd->sInwardCurrentBC.setZero();
  mpqd->ggqd->eInwardCurrentBC.setZero();

  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR; iR++)
  {
    // Set biases
    nBias = checkForZeroAxialFlux(0,iR);
    nInwardBias = checkForZeroInwardFluxNorthBC(iR);
    sBias = checkForZeroAxialFlux(mesh->nZ,iR);
    sInwardBias = checkForZeroInwardFluxSouthBC(iR);

    // Reset accumulators
    nFluxAccum = 0.0;
    sFluxAccum = 0.0;

    nInwardFluxAccum = 0.0;
    sInwardFluxAccum = 0.0;

    // Loop over neutron energy groups
    for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
    {

      // Get flux in these cells and energy group
      nFlux = mgqd->SGQDs[iEnergyGroup]->sFluxZ(0,iR) + nBias; 
      sFlux = mgqd->SGQDs[iEnergyGroup]->sFluxZ(mesh->nZ,iR) + sBias; 

      nInwardFlux = mgqd->SGQDs[iEnergyGroup]->nInwardFluxBC(iR) + nInwardBias; 
      sInwardFlux = mgqd->SGQDs[iEnergyGroup]->sInwardFluxBC(iR) + sInwardBias; 

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

      nInwardFluxAccum += nInwardFlux;
      sInwardFluxAccum += sInwardFlux;

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
      /nInwardFluxAccum; 
    mpqd->ggqd->sOutwardCurrToFluxRatioInwardWeightedBC(iR)\
      = mpqd->ggqd->sOutwardCurrToFluxRatioInwardWeightedBC(iR)\
      /sInwardFluxAccum;

    // Set flux boundaries
    mpqd->ggqd->nFluxBC(iR) = nFluxAccum - mats->nGroups*nBias;
    mpqd->ggqd->sFluxBC(iR) = sFluxAccum - mats->nGroups*sBias;

    mpqd->ggqd->nInwardFluxBC(iR) = nInwardFluxAccum - mats->nGroups*nInwardBias; 
    mpqd->ggqd->sInwardFluxBC(iR) = sInwardFluxAccum - mats->nGroups*sInwardBias; 

  }

  // Loop over spatial mesh
  for (int iZ = 0; iZ < mesh->nZ; iZ++)
  {

    // Set biases
    eBias = checkForZeroRadialFlux(iZ,mesh->nR);
    eInwardBias = checkForZeroInwardFluxEastBC(iZ);

    // Reset accumulators
    eFluxAccum = 0.0;

    eInwardFluxAccum = 0.0;

    // Loop over neutron energy groups
    for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
    {

      // Get flux in these cells and energy group
      eFlux = mgqd->SGQDs[iEnergyGroup]->sFluxR(iZ,mesh->nR) + eBias; 

      eInwardFlux = mgqd->SGQDs[iEnergyGroup]->eInwardFluxBC(iZ) + eInwardBias; 

      eInwardCurrent = mgqd->SGQDs[iEnergyGroup]->eInwardCurrentBC(iZ); 

      // Get outward current to flux ratio for these cells
      eRatio = mgqd->SGQDs[iEnergyGroup]->eOutwardCurrToFluxRatioBC(iZ); 

      // Flux weighted quasidiffusion coefficient
      mpqd->ggqd->eOutwardCurrToFluxRatioBC(iZ) += eRatio*eFlux;

      mpqd->ggqd->eOutwardCurrToFluxRatioInwardWeightedBC(iZ)\
        += eRatio*eInwardFlux;

      // Accumulate fluxes and currents 
      eFluxAccum = eFluxAccum + eFlux;

      eInwardFluxAccum += eInwardFlux;

      mpqd->ggqd->eInwardCurrentBC(iZ) += eInwardCurrent;
    }

    // Divide data by accumulated values
    mpqd->ggqd->eOutwardCurrToFluxRatioBC(iZ)\
      = mpqd->ggqd->eOutwardCurrToFluxRatioBC(iZ)/eFluxAccum; 

    mpqd->ggqd->eOutwardCurrToFluxRatioInwardWeightedBC(iZ)\
      = mpqd->ggqd->eOutwardCurrToFluxRatioInwardWeightedBC(iZ)\
      /eInwardFluxAccum;

    // Set flux boundaries
    mpqd->ggqd->eFluxBC(iZ) = eFluxAccum - mats->nGroups*eBias;

    mpqd->ggqd->eInwardFluxBC(iZ) = eInwardFluxAccum\
                                    - mats->nGroups*eInwardBias;
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
  double mySigT,myNeutV,eps;

  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ+1; iZ++)
    {

      // Check whether a bias is needed to prevent division by zero during
      // group collapse. eps = 0, if not. 
      eps = checkForZeroAxialCurrent(iZ,iR);           

      // Reset accumulator
      zCurrentAccum = 0.0;

      // Loop over neutron energy groups
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        zCurrent = abs(mgqd->SGQDs[iEnergyGroup]->currentZ(iZ,iR) + eps);

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
  double mySigT,myNeutV,eps;

  // Loop over spatial mesh
  for (int iR = 0; iR < mesh->nR+1; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ; iZ++)
    {

      // Check whether a bias is needed to prevent division by zero during
      // group collapse. eps = 0, if not. 
      eps = checkForZeroRadialCurrent(iZ,iR);           

      // Reset accumulators
      rCurrentAccum = 0.0;

      // Loop over neutron energy groups
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        // Get flux and currents in this cell and energy group
        rCurrent = abs(mgqd->SGQDs[iEnergyGroup]->currentR(iZ,iR) + eps); 

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
      //timeDerivative = (presentSum-pastSum)/mesh->dt;
      // The expression below is consistent with that used by Tamang, although
      // it seems the continuous expression would suggest the form commented 
      // above
      timeDerivative = (presentSum)/mesh->dt; 

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

        // Axial current weighted zeta1
        pastSum = pastSum + (1.0/myNeutV-1.0/myZNeutVPast)*zCurrentPast;
        presentSum = presentSum + (1.0/myNeutV-1.0/myZNeutV)*zCurrent;

        // Axial current weighted zeta2
        sigTDiff += (mySigT-mySigTR)*zCurrent;       

        // Accumulate absolute radial current
        fluxAccum = fluxAccum + flux;

      }

      // Approximate derivative
      //timeDerivative = (presentSum-pastSum)/mesh->dt; 
      // The expression below is consistent with that used by Tamang, although
      // it seems the continuous expression would suggest the form commented 
      // above
      timeDerivative = (presentSum)/mesh->dt; 

      // Divide by accumulated flux
      mats->oneGroupXS->zZeta1(iZ,iR) = timeDerivative/fluxAccum; 
      mats->oneGroupXS->zZeta2(iZ,iR) = sigTDiff/fluxAccum; 
      mats->oneGroupXS->zZeta(iZ,iR) = mats->oneGroupXS->zZeta1(iZ,iR)\
                                       + mats->oneGroupXS->zZeta2(iZ,iR); 
    }
  }
};
//==============================================================================

//==============================================================================
/// Check if flux is zero in every group and return a biasing factor 
/// to prevent division by zero   
///
/// @param [in] iZ axial cell index 
/// @param [in] iR radial cell index 
/// @param [out] eps bias factor 
double MGQDToMPQDCoupling::checkForZeroInwardFluxSouthBC(int iR)
{

  double groupFlux;
  bool aboveThreshold = true;  

  // Loop over energy groups 
  for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
  {
    // Check to see if current at location is greater than some eps
    groupFlux = mgqd->SGQDs[iEnergyGroup]->sInwardFluxBC(iR);
    aboveThreshold = aboveThreshold and (abs(groupFlux) > biasEps); 
  }

  if (aboveThreshold)
    return 0.0;
  else
    return biasEps;

};
//==============================================================================

//==============================================================================
/// Check if flux is zero in every group and return a biasing factor 
/// to prevent division by zero   
///
/// @param [in] iZ axial cell index 
/// @param [in] iR radial cell index 
/// @param [out] eps bias factor 
double MGQDToMPQDCoupling::checkForZeroInwardFluxNorthBC(int iR)
{

  double groupFlux;
  bool aboveThreshold = true;  

  // Loop over energy groups 
  for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
  {
    // Check to see if current at location is greater than some eps
    groupFlux = mgqd->SGQDs[iEnergyGroup]->nInwardFluxBC(iR);
    aboveThreshold = aboveThreshold and (abs(groupFlux) > biasEps); 
  }

  if (aboveThreshold)
    return 0.0;
  else
    return biasEps;

};
//==============================================================================

//==============================================================================
/// Check if flux is zero in every group and return a biasing factor 
/// to prevent division by zero   
///
/// @param [in] iZ axial cell index 
/// @param [in] iR radial cell index 
/// @param [out] eps bias factor 
double MGQDToMPQDCoupling::checkForZeroInwardFluxEastBC(int iZ)
{

  double groupFlux;
  bool aboveThreshold = true;  

  // Loop over energy groups 
  for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
  {
    // Check to see if current at location is greater than some eps
    groupFlux = mgqd->SGQDs[iEnergyGroup]->eInwardFluxBC(iZ);
    aboveThreshold = aboveThreshold and (abs(groupFlux) > biasEps); 
  }

  if (aboveThreshold)
    return 0.0;
  else
    return biasEps;

};
//==============================================================================

//==============================================================================
/// Check if flux is zero in every group and return a biasing factor 
/// to prevent division by zero   
///
/// @param [in] iZ axial cell index 
/// @param [in] iR radial cell index 
/// @param [out] eps bias factor 
double MGQDToMPQDCoupling::checkForZeroFlux(int iZ,int iR)
{

  double groupFlux;
  bool aboveThreshold = true;  

  // Loop over energy groups 
  for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
  {
    // Check to see if current at location is greater than some eps
    groupFlux = mgqd->SGQDs[iEnergyGroup]->sFlux(iZ,iR);
    aboveThreshold = aboveThreshold and (abs(groupFlux) > biasEps); 
  }

  if (aboveThreshold)
    return 0.0;
  else
    return biasEps;

};
//==============================================================================

//==============================================================================
/// Check if radial flux is zero in every group and return a biasing factor 
/// to prevent division by zero   
///
/// @param [in] iZ axial cell index 
/// @param [in] iR radial cell index 
/// @param [out] eps bias factor 
double MGQDToMPQDCoupling::checkForZeroRadialFlux(int iZ,int iR)
{

  double groupFlux;
  bool aboveThreshold = true;  

  // Loop over energy groups 
  for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
  {
    // Check to see if current at location is greater than some eps
    groupFlux = mgqd->SGQDs[iEnergyGroup]->sFluxR(iZ,iR);
    aboveThreshold = aboveThreshold and (abs(groupFlux) > biasEps); 
  }

  if (aboveThreshold)
    return 0.0;
  else
    return biasEps;

};
//==============================================================================

//==============================================================================
/// Check if axial flux is zero in every group and return a biasing factor 
/// to prevent division by zero   
///
/// @param [in] iZ axial cell index 
/// @param [in] iR radial cell index 
/// @param [out] eps bias factor 
double MGQDToMPQDCoupling::checkForZeroAxialFlux(int iZ,int iR)
{

  double groupFlux;
  bool aboveThreshold = true;  

  // Loop over energy groups 
  for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
  {
    // Check to see if current at location is greater than some eps
    groupFlux = mgqd->SGQDs[iEnergyGroup]->sFluxZ(iZ,iR);
    aboveThreshold = aboveThreshold and (abs(groupFlux) > biasEps); 
  }

  if (aboveThreshold)
    return 0.0;
  else
    return biasEps;

};
//==============================================================================

//==============================================================================
/// Check if axial current is zero in every group and return a biasing factor 
/// to prevent division by zero   
///
/// @param [in] iZ axial cell index 
/// @param [in] iR radial cell index 
/// @param [out] eps bias factor 
double MGQDToMPQDCoupling::checkForZeroAxialCurrent(int iZ,int iR)
{

  double groupCurrent;
  bool aboveThreshold = true;  

  // Loop over energy groups 
  for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
  {
    // Check to see if current at location is greater than some eps
    groupCurrent = mgqd->SGQDs[iEnergyGroup]->currentZ(iZ,iR);
    aboveThreshold = aboveThreshold and (abs(groupCurrent) > biasEps); 
  }

  if (aboveThreshold)
    return 0.0;
  else
    return biasEps;

};
//==============================================================================

//==============================================================================
/// Check if radial current is zero in every group and return a biasing factor 
/// to prevent division by zero   
///
/// @param [in] iZ axial cell index 
/// @param [in] iR radial cell index 
/// @param [out] eps bias factor 
double MGQDToMPQDCoupling::checkForZeroRadialCurrent(int iZ,int iR)
{

  double groupCurrent;
  bool aboveThreshold = true;  

  // Loop over energy groups 
  for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
  {
    // Check to see if current at location is greater than some eps
    groupCurrent = mgqd->SGQDs[iEnergyGroup]->currentR(iZ,iR);
    aboveThreshold = aboveThreshold and (abs(groupCurrent) > biasEps); 
  }

  if (aboveThreshold)
    return 0.0;
  else
    return biasEps;

};
//==============================================================================

//==============================================================================
/// Calculate residual between two vectors 
///
double MGQDToMPQDCoupling::calcResidual(Eigen::VectorXd vector1,\
    Eigen::VectorXd vector2)
{
  
  Eigen::VectorXd ones,residualVec,diff;
  double residual;
  double min,minScale = 1E-12;
  int fluxBeginIdx,fluxEndIdx,heatBeginIdx,heatEndIdx,dnpBeginIdx,dnpEndIdx; 

  // Set flux indices
  fluxBeginIdx = 0; 
  fluxEndIdx = mpqd->ggqd->nUnknowns;
  
  // Set temp indices
  heatBeginIdx = mpqd->heat->indexOffset;
  heatEndIdx = mpqd->heat->nUnknowns;
  
  // Set DNP indices
  dnpBeginIdx = mpqd->mgdnp->indexOffset;
  dnpEndIdx = vector1.size()-1;

  ones.setOnes(vector1.size());
 
  // Check for small values in flux 
  min = minScale*vector1(Eigen::seqN(fluxBeginIdx,fluxEndIdx)).mean();
  for (int index = 0; index < mpqd->ggqd->nUnknowns; index++)
  {
    if (abs(vector1(index)) < min and abs(vector2(index)) < min ) 
    {
      vector1(index) = 1.0;
      vector2(index) = 1.0;
    }
  }

  // Check for small values in temperature 
  min = minScale*vector1(Eigen::seqN(heatBeginIdx,heatEndIdx)).mean();
  for (int index = heatBeginIdx; index < dnpBeginIdx;\
      index++)
  {
    if (abs(vector1(index)) < min and abs(vector2(index)) < min ) 
    {
      vector1(index) = 1.0;
      vector2(index) = 1.0;
    }
  }

  // Check for small values in DNP concentrations 
  min = minScale*vector1(Eigen::seq(dnpBeginIdx,dnpEndIdx)).mean();
  for (int index = dnpBeginIdx; index < vector1.size(); index++)
  {
    if (abs(vector1(index)) < min and abs(vector2(index)) < min ) 
    {
      vector1(index) = 1.0;
      vector2(index) = 1.0;
    }
  }

  diff = (vector1-vector2);
  residualVec = diff.cwiseQuotient(vector1);
  residual = (1.0/residualVec.size())*residualVec.norm();
//  cout << "ResidualVec:" << endl;
//  cout << residualVec << endl;

  return residual;
  
};
//==============================================================================


