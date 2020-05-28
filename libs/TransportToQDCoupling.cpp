// File: TransportToQDCoupling.cpp
// Purpose: handle communication between transport and quasidiffusion objects
// Date: February 12, 2020

#include "TransportToQDCoupling.h"
#include "SimpleCornerBalance.h"

using namespace std;

//==============================================================================
/// TransportToQDCoupling class object constructor
///
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myInput YAML input object for the simulation
/// @param [in] myMGT Multigroup transport object for the simulation
/// @param [in] myMGQD Multigroup quasidiffusion object for the simulation
TransportToQDCoupling::TransportToQDCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput,\
  MultiGroupTransport * myMGT,\
  MultiGroupQD * myMGQD)
{
  // Assign pointers for materials, mesh, input objects, multigroup transport
  // and quasidiffusion objects
  materials = myMaterials;
  mesh = myMesh;
  input = myInput;
  MGT = myMGT;
  MGQD = myMGQD;
  
  // Check for optional parameters
  checkOptionalParams(); 

};
//==============================================================================

//==============================================================================
/// Calculate Eddington factors using angular fluxes from transport objects
/// @param [out] allConverged boolean indicating if the residual on the 
/// Eddington factors passes the convergence criteria
bool TransportToQDCoupling::calcEddingtonFactors()
{
  int rows = MGT->SGTs[0]->sFlux.rows();
  int cols = MGT->SGTs[0]->sFlux.cols();
  int angIdx,xiIdx=0,muIdx=1,etaIdx=2,weightIdx = 3;  
  double angFlux,mu,xi,weight,EzzCoef,ErrCoef,ErzCoef;
  double numeratorEzz,numeratorErr,numeratorErz,denominator;
  double residualZz,residualRr,residualRz;
  bool allConverged=true;
  
  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    // store past eddington factors
    MGQD->SGQDs[iGroup]->EzzPrev = MGQD->SGQDs[iGroup]->Ezz;
    MGQD->SGQDs[iGroup]->ErrPrev = MGQD->SGQDs[iGroup]->Err;
    MGQD->SGQDs[iGroup]->ErzPrev = MGQD->SGQDs[iGroup]->Erz;

    for (int iR = 0; iR < cols; iR++)
    {
      for (int iZ = 0; iZ < rows; iZ++)
      {

        // reset accumulators
        numeratorEzz = 0.0;
        numeratorErr = 0.0;
        numeratorErz = 0.0;
        denominator = 0.0;

        // loop over quadrature
        for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi)
        {
          xi = mesh->quadrature[iXi].quad[0][xiIdx];
          for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu)
          {
            angIdx=mesh->quadrature[iXi].ordIdx[iMu];
            mu = mesh->quadrature[iXi].quad[iMu][muIdx]; 
            weight = mesh->quadrature[iXi].quad[iMu][weightIdx]; 
            angFlux = MGT->SGTs[iGroup]->aFlux(iZ,iR,angIdx);

            EzzCoef = xi*xi;
            ErrCoef = mu*mu;
            ErzCoef = mu*xi;

            numeratorEzz = numeratorEzz + EzzCoef*angFlux*weight;
            numeratorErr = numeratorErr + ErrCoef*angFlux*weight;
            numeratorErz = numeratorErz + ErzCoef*angFlux*weight;

            denominator = denominator + angFlux*weight;
          
          } //iMu
        } //iXi 
 
        MGQD->SGQDs[iGroup]->Ezz(iZ,iR) = numeratorEzz/denominator;
        MGQD->SGQDs[iGroup]->Err(iZ,iR) = numeratorErr/denominator;
        MGQD->SGQDs[iGroup]->Erz(iZ,iR) = numeratorErz/denominator;

      } //iZ
    } //iR
   
    // measure the residual for each Eddington factors 
    residualZz = ((MGQD->SGQDs[iGroup]->Ezz - MGQD->SGQDs[iGroup]->EzzPrev)\
      .cwiseQuotient(MGQD->SGQDs[iGroup]->Ezz)).norm();
    residualRr = ((MGQD->SGQDs[iGroup]->Err - MGQD->SGQDs[iGroup]->ErrPrev)\
      .cwiseQuotient(MGQD->SGQDs[iGroup]->Err)).norm();
    residualRz = ((MGQD->SGQDs[iGroup]->Erz - MGQD->SGQDs[iGroup]->ErzPrev)\
      .cwiseQuotient(MGQD->SGQDs[iGroup]->Erz)).norm();
    cout << "eddington residual" << endl; 
    cout << residualZz << endl;
    cout << residualRr << endl;
    cout << residualRz << endl;
    
    if (residualZz < epsEddington and residualRr < epsEddington and 
      residualRz < epsEddington)
      allConverged = allConverged and true;
    else
      allConverged = allConverged and false;
  } //iGroup

  return allConverged;
}
//==============================================================================

//==============================================================================
/// Calculate a number a parameters used for forming the boundary conditions of
/// low order problem 
void TransportToQDCoupling::calcBCs()
{
  int rows = MGT->SGTs[0]->sFlux.rows();
  int cols = MGT->SGTs[0]->sFlux.cols();
  int eIdx = cols - 1,sIdx = rows - 1;
  int angIdx,xiIdx=0,muIdx=1,etaIdx=2,weightIdx = 3;  
  double angFlux,angFluxN,angFluxS,mu,xi,weight;
  double inwardJrE,inwardJzN,inwardJzS;
  double inwardFluxE,inwardFluxN,inwardFluxS;
  double outwardJrE,outwardJzN,outwardJzS;
  double outwardFluxE,outwardFluxN,outwardFluxS;
  double localScalarFluxE, localScalarFluxN, localScalarFluxS;
  
  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
      for (int iZ = 0; iZ < rows; iZ++)
      {
        // reset accumulators
        inwardJrE = 0.0;
        inwardFluxE = 0.0;
        outwardJrE = 0.0;
        outwardFluxE = 0.0;
        localScalarFluxE = 0.0;

        // loop over quadrature
        for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi)
        {
          xi = mesh->quadrature[iXi].quad[0][xiIdx];
          for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu)
          {
            angIdx=mesh->quadrature[iXi].ordIdx[iMu];
            mu = mesh->quadrature[iXi].quad[iMu][muIdx]; 
            weight = mesh->quadrature[iXi].quad[iMu][weightIdx];
 
            angFlux = MGT->SGTs[iGroup]->aFlux(iZ,eIdx,angIdx);
            
            localScalarFluxE += angFlux*weight;

            // only accumulate inward facing angular fluxes on the the 
            // outside radial (east) boundary 
            if (mu < 0) 
            {
              inwardJrE = inwardJrE + mu*angFlux*weight;
              inwardFluxE = inwardFluxE + angFlux*weight;
            } else
            {
              outwardJrE = outwardJrE + mu*angFlux*weight;
              outwardFluxE = outwardFluxE + angFlux*weight;
            }
          } //iMu
        } //iXi 

        // set inward current in SGQD object 
        MGQD->SGQDs[iGroup]->eInwardCurrentBC(iZ) = inwardJrE;

        // set inward flux in SGQD object 
        MGQD->SGQDs[iGroup]->eInwardFluxBC(iZ) = inwardFluxE;

        // set outward current to flux ratio in SGQD object 
        MGQD->SGQDs[iGroup]->eOutwardCurrToFluxRatioBC(iZ)\
          = outwardJrE/outwardFluxE;

        // set eastern flux bc
        MGQD->SGQDs[iGroup]->eFluxBC(iZ) = localScalarFluxE; 
        MGQD->SGQDs[iGroup]->eCurrentRBC(iZ) = outwardJrE+inwardJrE;

        // set absolute current
        MGQD->SGQDs[iGroup]->eAbsCurrentBC(iZ) = outwardJrE - inwardJrE;
      
        } //iZ
      for (int iR = 0; iR < cols; iR++)
      {
        // reset accumulators
        inwardJzN = 0.0;
        inwardJzS = 0.0;
        inwardFluxN = 0.0;
        inwardFluxS = 0.0;
        outwardJzN = 0.0;
        outwardJzS = 0.0;
        outwardFluxN = 0.0;
        outwardFluxS = 0.0;
        localScalarFluxN = 0.0;
        localScalarFluxS = 0.0;

        // loop over quadrature
        for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi)
        {
          xi = mesh->quadrature[iXi].quad[0][xiIdx];
          for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu)
          {
            angIdx=mesh->quadrature[iXi].ordIdx[iMu];
            mu = mesh->quadrature[iXi].quad[iMu][muIdx]; 
            weight = mesh->quadrature[iXi].quad[iMu][weightIdx];
 
            angFluxN = MGT->SGTs[iGroup]->aFlux(0,iR,angIdx);
            angFluxS = MGT->SGTs[iGroup]->aFlux(sIdx,iR,angIdx);
  
            localScalarFluxN += angFluxN*weight;
            localScalarFluxS += angFluxS*weight;

            // accumulate outward and inward angular fluxes in separate
            // variables on the north face
            if (xi > 0) 
            {
              inwardJzN = inwardJzN + xi*angFluxN*weight;
              inwardFluxN = inwardFluxN + angFluxN*weight;
            } else
            {
              outwardJzN = outwardJzN + xi*angFluxN*weight;
              outwardFluxN = outwardFluxN + angFluxN*weight;
            }
            // accumulate outward and inward angular fluxes in separate
            // variables on the south face
            if (xi < 0)
            {
              inwardJzS = inwardJzS + xi*angFluxS*weight;
              inwardFluxS = inwardFluxS + angFluxS*weight;
            } else
            {
              outwardJzS = outwardJzS + xi*angFluxS*weight;
              outwardFluxS = outwardFluxS + angFluxS*weight;
            }
          
          } //iMu
        } //iXi 

        // set inward current in SGQD object 
        MGQD->SGQDs[iGroup]->nInwardCurrentBC(iR) = inwardJzN;
        MGQD->SGQDs[iGroup]->sInwardCurrentBC(iR) = inwardJzS;
        
        // set inward flux in SGQD object 
        MGQD->SGQDs[iGroup]->nInwardFluxBC(iR) = inwardFluxN;
        MGQD->SGQDs[iGroup]->sInwardFluxBC(iR) = inwardFluxS;
  
        // set outward current to flux ratio in SGQD object
        MGQD->SGQDs[iGroup]->nOutwardCurrToFluxRatioBC(iR)\
          = outwardJzN/outwardFluxN;
        MGQD->SGQDs[iGroup]->sOutwardCurrToFluxRatioBC(iR)\
          = outwardJzS/outwardFluxS;
        
        // set flux BCs
        MGQD->SGQDs[iGroup]->nFluxBC(iR) = localScalarFluxN;
        MGQD->SGQDs[iGroup]->sFluxBC(iR) = localScalarFluxS;
        MGQD->SGQDs[iGroup]->nCurrentZBC(iR) = outwardJzN+inwardJzN;
        MGQD->SGQDs[iGroup]->sCurrentZBC(iR) = outwardJzS+inwardJzS;
        
        // set absolute currents
        MGQD->SGQDs[iGroup]->nAbsCurrentBC(iR) = outwardJzN - inwardJzN;
        MGQD->SGQDs[iGroup]->sAbsCurrentBC(iR) = outwardJzS - inwardJzS;

      } //iR 
  
  } //iGroup

}
//==============================================================================

//==============================================================================
/// Use a multigroup transport solve to form the BCs and Eddington factors,
/// which are then used in a multigroup quasidiffusion solve. The solution of 
/// the MGQD system is used to update the sources in the transport problem. This
/// process repeats until the eddington factors, sources and alphas are 
/// converged.
void TransportToQDCoupling::solveTransportWithQDAcceleration()
{

  bool alphaConverged=false,eddingtonConverged=false,sourcesConverged=false;

  MGQD->setInitialCondition();
  
  // loop over time steps
  for (int iTime = 0; iTime < mesh->dts.size(); iTime++)
  {
 
    MGQD->buildLinearSystem();
    MGQD->solveLinearSystem();
    updateTransportFluxes();
    MGT->calcSources("fs");

    // iterate until eddington and alpha are converged. Max iterations set
    // to 1000 for now. ToDo: make maxIters and input.
    for (int iSteps = 0; iSteps < 1000; iSteps++)
    {
      MGT->solveStartAngles();
      MGT->solveSCBs();
      MGT->calcFluxes();
      eddingtonConverged = calcEddingtonFactors();
      calcBCs();
      MGQD->buildLinearSystem();
      MGQD->solveLinearSystem();
      updateTransportFluxes();
      sourcesConverged = MGT->calcSources("fs");
      if (sourcesConverged) alphaConverged = MGT->calcAlphas("print");
      if (alphaConverged and eddingtonConverged and sourcesConverged) break;
    }
  
    MGQD->writeFluxes();
    MGT->calcFluxes();
    MGT->writeFluxes();
    MGQD->QDSolve->xPast = MGQD->QDSolve->x;
    MGQD->buildBackCalcSystem();
    MGQD->backCalculateCurrent();
    updateTransportPrevFluxes();
  }
}
//==============================================================================

//==============================================================================
/// Update the transport fluxes with those in the quasidiffusion objects
void TransportToQDCoupling::updateTransportFluxes()
{
  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    MGT->SGTs[iGroup]->sFlux = MGQD->SGQDs[iGroup]->sFlux;
  }
}
//=============================================================================

//==============================================================================
/// Update transport fluxes at the previous time step with the fluxes currently
/// on the quasidiffusion objects
void TransportToQDCoupling::updateTransportPrevFluxes()
{
  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    MGT->SGTs[iGroup]->sFluxPrev = MGQD->SGQDs[iGroup]->sFlux;
  }
}
//=============================================================================

//==============================================================================
/// Check for optional input parameters of relevance to this object
void TransportToQDCoupling::checkOptionalParams()
{
  if ((*input)["parameters"]["epsEddington"])
  {
    epsEddington=(*input)["parameters"]["epsEddington"].as<double>();
  }
}
//==============================================================================
