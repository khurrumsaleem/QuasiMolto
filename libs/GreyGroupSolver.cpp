// File: GreyGroupSolver.cpp
// Purpose: Solve RZ quasidiffusion equations  
// Date: February 05, 2020

#include "GreyGroupSolver.h"
#include "GreyGroupQD.h"

using namespace std; 

//==============================================================================
/// GreyGroupSolver object constructor
///
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myInput YAML input object for the simulation
/// @param [in] myMPQD multiphysics coupled QD object for the simulation
GreyGroupSolver::GreyGroupSolver(GreyGroupQD * myGGQD,\
    Mesh * myMesh,\
    Materials * myMaterials,\
    YAML::Node * myInput)	      
{
  // Point to variables for mesh and input file
  GGQD = myGGQD;
  mesh = myMesh;
  input = myInput;
  materials = myMaterials;

  // calculate number of unknowns  
  nR = mesh->nR;
  nZ = mesh->nZ;
  nUnknowns = 3*(nZ*nR) + nZ + nR;
  nCurrentUnknowns = 2*(nZ*nR) + nZ + nR;

  // initialize size of linear system
  C.resize(nCurrentUnknowns,nUnknowns);
  C.reserve(4*nCurrentUnknowns);
  xFlux.setZero(nUnknowns);
  currPast.setZero(nCurrentUnknowns);
  d.setZero(nCurrentUnknowns);

  // Set number of unknowns if GreyGroupQD object
  GGQD->nUnknowns = nUnknowns;

  checkOptionalParams();
};

//==============================================================================

//==============================================================================
/// Form a portion of the linear system  

void GreyGroupSolver::formLinearSystem()	      
{

  int iEq = GGQD->indexOffset;
  Atemp.resize(nUnknowns,A->cols());
  Atemp.reserve(10*nUnknowns);

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {

      // apply zeroth moment equation
      assertZerothMoment(iR,iZ,iEq);
      iEq = iEq + 1;

      // south face
      if (iZ == mesh->dzsCorner.size()-1)
      {
        // if on the boundary, assert boundary conditions
        assertSBC(iR,iZ,iEq);
        iEq = iEq + 1;
      } else
      {
        // otherwise assert first moment balance on south face
        applyAxialBoundary(iR,iZ,iEq);
        iEq = iEq + 1;
      }

      // east face
      if (iR == mesh->drsCorner.size()-1)
      {
        // if on the boundary, assert boundary conditions
        assertEBC(iR,iZ,iEq);
        iEq = iEq + 1;
      } else
      {
        // otherwise assert first moment balance on north face
        applyRadialBoundary(iR,iZ,iEq);
        iEq = iEq + 1;
      }

      // north face
      if (iZ == 0)
      {
        // if on the boundary, assert boundary conditions
        assertNBC(iR,iZ,iEq);
        iEq = iEq + 1;
      } 

      // west face
      if (iR == 0)
      {
        // if on the boundary, assert boundary conditions
        assertWBC(iR,iZ,iEq);
        iEq = iEq + 1;
      } 

    }
  }

  A->middleRows(GGQD->indexOffset,nUnknowns) = Atemp; 
};

//==============================================================================

//==============================================================================
/// Form a portion of the linear system  

void GreyGroupSolver::formSteadyStateLinearSystem()	      
{

  int iEq = GGQD->indexOffset;
  Atemp.resize(nUnknowns,A->cols());
  Atemp.reserve(10*nUnknowns);

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {

      // apply zeroth moment equation
      assertSteadyStateZerothMoment(iR,iZ,iEq);
      iEq = iEq + 1;

      // south face
      if (iZ == mesh->dzsCorner.size()-1)
      {
        // if on the boundary, assert boundary conditions
        assertSteadyStateSBC(iR,iZ,iEq);
        iEq = iEq + 1;
      } else
      {
        // otherwise assert first moment balance on south face
        applySteadyStateAxialBoundary(iR,iZ,iEq);
        iEq = iEq + 1;
      }

      // east face
      if (iR == mesh->drsCorner.size()-1)
      {
        // if on the boundary, assert boundary conditions
        assertSteadyStateEBC(iR,iZ,iEq);
        iEq = iEq + 1;
      } else
      {
        // otherwise assert first moment balance on north face
        applySteadyStateRadialBoundary(iR,iZ,iEq);
        iEq = iEq + 1;
      }

      // north face
      if (iZ == 0)
      {
        // if on the boundary, assert boundary conditions
        assertSteadyStateNBC(iR,iZ,iEq);
        iEq = iEq + 1;
      } 

      // west face
      if (iR == 0)
      {
        // if on the boundary, assert boundary conditions
        assertSteadyStateWBC(iR,iZ,iEq);
        iEq = iEq + 1;
      } 

    }
  }

  A->middleRows(GGQD->indexOffset,nUnknowns) = Atemp; 
};

//==============================================================================

//==============================================================================
/// Compute currents from flux values in x
///
void GreyGroupSolver::backCalculateCurrent()
{
  currPast = d + C*(xFlux); 
}
//==============================================================================

//==============================================================================
/// Assert the transient zeroth moment equation for cell (iR,iZ)
///
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertZerothMoment(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double vPast = materials->oneGroupXS->neutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double groupSourceCoeff;

  indices = getIndices(iR,iZ);

  // Scatter and fission source term
  groupSourceCoeff = calcScatterAndFissionCoeff(iR,iZ);
  Atemp.insert(iEq,indices[iCF]) = -geoParams[iCF] * groupSourceCoeff;

  // DNP source term
  GGQD->mpqd->dnpSource(iZ,iR,iEq,-geoParams[iCF], &Atemp);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  Atemp.coeffRef(iEq,indices[iCF]) += geoParams[iCF] * ((1/(v*deltaT)) + sigT);

  westCurrent(-geoParams[iWF],iR,iZ,iEq);

  eastCurrent(geoParams[iEF],iR,iZ,iEq);

  northCurrent(-geoParams[iNF],iR,iZ,iEq);

  southCurrent(geoParams[iSF],iR,iZ,iEq);

  // formulate RHS entry
  (*b)(iEq) = (*b)(iEq) + geoParams[iCF]*\
              ( ((*xPast)(indices[iCF])/(vPast*deltaT)) + GGQD->q(iZ,iR));
};
//==============================================================================


//==============================================================================
/// Apply transient radial boundary for cell (iR,iZ)
///
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::applyRadialBoundary(int iR,int iZ,int iEq)
{
  eastCurrent(1,iR,iZ,iEq);
  westCurrent(-1,iR+1,iZ,iEq);
}
//==============================================================================

//==============================================================================
/// Apply transient axial boundary for cell (iR,iZ)
///
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::applyAxialBoundary(int iR,int iZ,int iEq)
{
  northCurrent(1,iR,iZ+1,iEq);
  southCurrent(-1,iR,iZ,iEq);
}
//==============================================================================

//==============================================================================
/// Assert the steady state zeroth moment equation for cell (iR,iZ)
///
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateZerothMoment(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double vPast = materials->oneGroupXS->neutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double scatterCoeff, fissionCoeff, keff, cellFlux;

  indices = getIndices(iR,iZ);

  // Scattering source term (implicit)
  scatterCoeff = materials->oneGroupXS->sigS(iZ,iR);
  Atemp.insert(iEq,indices[iCF]) = -geoParams[iCF] * scatterCoeff;

  // Fission source term (explicit for power iteration)
  fissionCoeff = materials->oneGroupXS->qdFluxCoeff(iZ,iR);
  keff = materials->oneGroupXS->keff;
  cellFlux = GGQD->sFlux(iZ,iR);
  (*b)(iEq) = (*b)(iEq) + geoParams[iCF]*\
              ( fissionCoeff*cellFlux/keff + GGQD->q(iZ,iR));

  // DNP source term
  GGQD->mpqd->dnpSource(iZ,iR,iEq,-geoParams[iCF], &Atemp);

  // populate entries representing streaming and reaction terms
  Atemp.coeffRef(iEq,indices[iCF]) += geoParams[iCF] * (sigT);

  steadyStateWestCurrent(-geoParams[iWF],iR,iZ,iEq);

  steadyStateEastCurrent(geoParams[iEF],iR,iZ,iEq);

  steadyStateNorthCurrent(-geoParams[iNF],iR,iZ,iEq);

  steadyStateSouthCurrent(geoParams[iSF],iR,iZ,iEq);

};
//==============================================================================

//==============================================================================
/// Apply steady state radial boundary for cell (iR,iZ)
///
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::applySteadyStateRadialBoundary(int iR,int iZ,int iEq)
{
  steadyStateEastCurrent(1,iR,iZ,iEq);
  steadyStateWestCurrent(-1,iR+1,iZ,iEq);
}
//==============================================================================

//==============================================================================
/// Apply steady state axial boundary for cell (iR,iZ)
///
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::applySteadyStateAxialBoundary(int iR,int iZ,int iEq)
{
  steadyStateNorthCurrent(1,iR,iZ+1,iEq);
  steadyStateSouthCurrent(-1,iR,iZ,iEq);
}
//==============================================================================

//==============================================================================
/// Enforce coefficients for transient current on south face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::southCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->zNeutV(iZ+1,iR);
  double vPast = materials->oneGroupXS->zNeutVPast(iZ+1,iR);
  double sigT = materials->oneGroupXS->zSigTR(iZ+1,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,zetaL,\
    mgqdCurrent,mgqdNeutV;  
  double EzzC,EzzS,ErzW,ErzE;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zUp-zAvg;

  // get local Eddington factors
  EzzC = GGQD->Ezz(iZ,iR);
  EzzS = GGQD->EzzAxial(iZ+1,iR); 
  ErzW = GGQD->ErzRadial(iZ,iR);
  ErzE = GGQD->ErzRadial(iZ,iR+1);
  zetaL = materials->oneGroupXS->zZeta(iZ+1,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = coeff/((1/(v*deltaT))+sigT); 

  Atemp.coeffRef(iEq,indices[iSF]) -= coeff*EzzS/deltaZ;

  Atemp.coeffRef(iEq,indices[iCF]) += coeff*EzzC/deltaZ;

  Atemp.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzW/(rAvg*deltaR));

  Atemp.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzE/(rAvg*deltaR);

  // Enforce zeta coefficient
  Atemp.coeffRef(iEq,indices[iSF]) -= coeff*zetaL;

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      mgqdCurrent = GGQD->mgqd->QDSolve->currPast(indices[iSC]); 
      mgqdNeutV = materials->zNeutVel(iZ+1,iR,iGroup); 
      (*b)(iEq) -= (coeff/deltaT)*(mgqdCurrent/mgqdNeutV);
    }
  }
  else
  {
    (*b)(iEq) = (*b)(iEq) - coeff*(currPast(indices[iSC])/(vPast*deltaT));
  }

};
//==============================================================================

//==============================================================================
/// Enforce coefficients for transient current on north face 
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::northCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->zNeutV(iZ,iR);
  double vPast = materials->oneGroupXS->zNeutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->zSigTR(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,zetaL,\
    mgqdCurrent,mgqdNeutV;  
  double EzzC,EzzN,ErzW,ErzE;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zAvg-zDown;

  // get local Eddington factors
  EzzC = GGQD->Ezz(iZ,iR);
  EzzN = GGQD->EzzAxial(iZ,iR); 
  ErzW = GGQD->ErzRadial(iZ,iR);
  ErzE = GGQD->ErzRadial(iZ,iR+1);
  zetaL = materials->oneGroupXS->zZeta(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = coeff/((1/(v*deltaT))+sigT); 

  Atemp.coeffRef(iEq,indices[iNF]) += coeff*EzzN/deltaZ;

  Atemp.coeffRef(iEq,indices[iCF]) -= coeff*EzzC/deltaZ;

  Atemp.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzW/(rAvg*deltaR));

  Atemp.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzE/(rAvg*deltaR);

  // Enforce zeta coefficient
  Atemp.coeffRef(iEq,indices[iNF]) -= coeff*zetaL;

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      mgqdCurrent = GGQD->mgqd->QDSolve->currPast(indices[iNC]); 
      mgqdNeutV = materials->zNeutVel(iZ,iR,iGroup); 
      (*b)(iEq) -= (coeff/deltaT)*(mgqdCurrent/mgqdNeutV);
    }
  }
  else
  {
    (*b)(iEq) = (*b)(iEq) - coeff*(currPast(indices[iNC])/(vPast*deltaT));
  }
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for transient current on west face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::westCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->rNeutV(iZ,iR);
  double vPast = materials->oneGroupXS->rNeutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->rSigTR(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,zetaL;
  double hCent,hDown,mgqdCurrent,mgqdNeutV;
  double ErzN,ErzS,ErrC,ErrW;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rAvg-rDown; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg,iWF);
  hDown = calcIntegratingFactor(iR,iZ,rDown,iWF);

  // get local Eddington factors
  ErrC = GGQD->Err(iZ,iR);
  ErrW = GGQD->ErrRadial(iZ,iR); 
  ErzN = GGQD->ErzAxial(iZ,iR);
  ErzS = GGQD->ErzAxial(iZ+1,iR);
  zetaL = materials->oneGroupXS->rZeta(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = coeff/((1/(v*deltaT))+sigT); 

  Atemp.coeffRef(iEq,indices[iSF]) -= coeff*ErzS/deltaZ;

  Atemp.coeffRef(iEq,indices[iNF]) += coeff*ErzN/deltaZ;

  Atemp.coeffRef(iEq,indices[iCF]) -= coeff*hCent*ErrC/(hDown*deltaR);

  Atemp.coeffRef(iEq,indices[iWF]) += coeff*hDown*ErrW/(hDown*deltaR);

  // Enforce zeta coefficient
  Atemp.coeffRef(iEq,indices[iWF]) -= coeff*zetaL;

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      mgqdCurrent = GGQD->mgqd->QDSolve->currPast(indices[iWC]); 
      mgqdNeutV = materials->rNeutVel(iZ,iR,iGroup); 
      (*b)(iEq) -= (coeff/deltaT)*(mgqdCurrent/mgqdNeutV);
    }
  }
  else
  {
    (*b)(iEq) = (*b)(iEq) - coeff*(currPast(indices[iWC])/(v*deltaT));
  }

};
//==============================================================================

//==============================================================================
/// Enforce coefficients for transient current on east face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::eastCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->rNeutV(iZ,iR+1);
  double vPast = materials->oneGroupXS->rNeutVPast(iZ,iR+1);
  double sigT = materials->oneGroupXS->rSigTR(iZ,iR+1);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,zetaL;  
  double hCent,hUp,mgqdCurrent,mgqdNeutV;
  double ErzN,ErzS,ErrC,ErrE;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rAvg; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg,iEF);
  hUp = calcIntegratingFactor(iR,iZ,rUp,iEF);

  // get local Eddington factors
  ErrC = GGQD->Err(iZ,iR);
  ErrE = GGQD->ErrRadial(iZ,iR+1); 
  ErzN = GGQD->ErzAxial(iZ,iR);
  ErzS = GGQD->ErzAxial(iZ+1,iR);
  zetaL = materials->oneGroupXS->rZeta(iZ,iR+1);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = coeff/((1/(v*deltaT))+sigT); 

  Atemp.coeffRef(iEq,indices[iSF]) -= coeff*ErzS/deltaZ;

  Atemp.coeffRef(iEq,indices[iNF]) += coeff*ErzN/deltaZ;

  Atemp.coeffRef(iEq,indices[iCF]) += coeff*hCent*ErrC/(hUp*deltaR);

  Atemp.coeffRef(iEq,indices[iEF]) -= coeff*hUp*ErrE/(hUp*deltaR);

  // Enforce zeta coefficient
  Atemp.coeffRef(iEq,indices[iEF]) -= coeff*zetaL;

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      mgqdCurrent = GGQD->mgqd->QDSolve->currPast(indices[iEC]); 
      mgqdNeutV = materials->rNeutVel(iZ,iR+1,iGroup); 
      (*b)(iEq) -= (coeff/deltaT)*(mgqdCurrent/mgqdNeutV);
    }
  }
  else
  {
    (*b)(iEq) = (*b)(iEq) - coeff*(currPast(indices[iEC])/(vPast*deltaT));  
  }

};
//==============================================================================

//==============================================================================
/// Enforce coefficients for steady state current on south face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::steadyStateSouthCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->zNeutV(iZ+1,iR);
  double vPast = materials->oneGroupXS->zNeutVPast(iZ+1,iR);
  double sigT = materials->oneGroupXS->zSigTR(iZ+1,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,zetaL,\
    mgqdCurrent,mgqdNeutV;  
  double EzzC,EzzS,ErzW,ErzE;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zUp-zAvg;

  // get local Eddington factors
  EzzC = GGQD->Ezz(iZ,iR);
  EzzS = GGQD->EzzAxial(iZ+1,iR); 
  ErzW = GGQD->ErzRadial(iZ,iR);
  ErzE = GGQD->ErzRadial(iZ,iR+1);
  zetaL = materials->oneGroupXS->zZeta2(iZ+1,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = coeff/(sigT); 

  Atemp.coeffRef(iEq,indices[iSF]) -= coeff*EzzS/deltaZ;

  Atemp.coeffRef(iEq,indices[iCF]) += coeff*EzzC/deltaZ;

  Atemp.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzW/(rAvg*deltaR));

  Atemp.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzE/(rAvg*deltaR);

  // Enforce zeta coefficient
  Atemp.coeffRef(iEq,indices[iSF]) -= coeff*zetaL;

};
//==============================================================================

//==============================================================================
/// Enforce coefficients for steady state current on north face 
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::steadyStateNorthCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->zNeutV(iZ,iR);
  double vPast = materials->oneGroupXS->zNeutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->zSigTR(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,zetaL,\
    mgqdCurrent,mgqdNeutV;  
  double EzzC,EzzN,ErzW,ErzE;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zAvg-zDown;

  // get local Eddington factors
  EzzC = GGQD->Ezz(iZ,iR);
  EzzN = GGQD->EzzAxial(iZ,iR); 
  ErzW = GGQD->ErzRadial(iZ,iR);
  ErzE = GGQD->ErzRadial(iZ,iR+1);
  zetaL = materials->oneGroupXS->zZeta2(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = coeff/(sigT); 

  Atemp.coeffRef(iEq,indices[iNF]) += coeff*EzzN/deltaZ;

  Atemp.coeffRef(iEq,indices[iCF]) -= coeff*EzzC/deltaZ;

  Atemp.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzW/(rAvg*deltaR));

  Atemp.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzE/(rAvg*deltaR);

  // Enforce zeta coefficient
  Atemp.coeffRef(iEq,indices[iNF]) -= coeff*zetaL;

};
//==============================================================================

//==============================================================================
/// Enforce coefficients for steady state current on west face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::steadyStateWestCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->rNeutV(iZ,iR);
  double vPast = materials->oneGroupXS->rNeutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->rSigTR(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,zetaL;
  double hCent,hDown,mgqdCurrent,mgqdNeutV;
  double ErzN,ErzS,ErrC,ErrW;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rAvg-rDown; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg,iWF);
  hDown = calcIntegratingFactor(iR,iZ,rDown,iWF);

  // get local Eddington factors
  ErrC = GGQD->Err(iZ,iR);
  ErrW = GGQD->ErrRadial(iZ,iR); 
  ErzN = GGQD->ErzAxial(iZ,iR);
  ErzS = GGQD->ErzAxial(iZ+1,iR);
  zetaL = materials->oneGroupXS->rZeta2(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = coeff/(sigT); 

  Atemp.coeffRef(iEq,indices[iSF]) -= coeff*ErzS/deltaZ;

  Atemp.coeffRef(iEq,indices[iNF]) += coeff*ErzN/deltaZ;

  Atemp.coeffRef(iEq,indices[iCF]) -= coeff*hCent*ErrC/(hDown*deltaR);

  Atemp.coeffRef(iEq,indices[iWF]) += coeff*hDown*ErrW/(hDown*deltaR);

  // Enforce zeta coefficient
  Atemp.coeffRef(iEq,indices[iWF]) -= coeff*zetaL;

};
//==============================================================================

//==============================================================================
/// Enforce coefficients for steady state current on east face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::steadyStateEastCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->rNeutV(iZ,iR+1);
  double vPast = materials->oneGroupXS->rNeutVPast(iZ,iR+1);
  double sigT = materials->oneGroupXS->rSigTR(iZ,iR+1);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,zetaL;  
  double hCent,hUp,mgqdCurrent,mgqdNeutV;
  double ErzN,ErzS,ErrC,ErrE;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rAvg; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg,iEF);
  hUp = calcIntegratingFactor(iR,iZ,rUp,iEF);

  // get local Eddington factors
  ErrC = GGQD->Err(iZ,iR);
  ErrE = GGQD->ErrRadial(iZ,iR+1); 
  ErzN = GGQD->ErzAxial(iZ,iR);
  ErzS = GGQD->ErzAxial(iZ+1,iR);
  zetaL = materials->oneGroupXS->rZeta2(iZ,iR+1);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = coeff/(sigT); 

  Atemp.coeffRef(iEq,indices[iSF]) -= coeff*ErzS/deltaZ;

  Atemp.coeffRef(iEq,indices[iNF]) += coeff*ErzN/deltaZ;

  Atemp.coeffRef(iEq,indices[iCF]) += coeff*hCent*ErrC/(hUp*deltaR);

  Atemp.coeffRef(iEq,indices[iEF]) -= coeff*hUp*ErrE/(hUp*deltaR);

  // Enforce zeta coefficient
  Atemp.coeffRef(iEq,indices[iEF]) -= coeff*zetaL;

};
//==============================================================================

//==============================================================================
/// Form a portion of the current back calc linear system that belongs to GGQD 
///
void GreyGroupSolver::formBackCalcSystem()	      
{
  int iEq = GGQD->indexOffset;

  // Reset linear system
  C.setZero();
  d.setZero();

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {

      // south face
      calcSouthCurrent(iR,iZ,iEq);
      iEq = iEq + 1;

      // east face
      calcEastCurrent(iR,iZ,iEq);
      iEq = iEq + 1;

      // north face
      if (iZ == 0)
      {
        calcNorthCurrent(iR,iZ,iEq);
        iEq = iEq + 1;
      } 

      // west face
      if (iR == 0)
      {
        // if on the boundary, assert boundary conditions
        calcWestCurrent(iR,iZ,iEq);
        iEq = iEq + 1;
      } 

    }
  }
};

//==============================================================================

//==============================================================================
/// Form a portion of the current back calc linear system that belongs to GGQD 
///
void GreyGroupSolver::formSteadyStateBackCalcSystem()	      
{
  int iEq = GGQD->indexOffset;

  // Reset linear system
  C.setZero();
  d.setZero();

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {

      // south face
      calcSteadyStateSouthCurrent(iR,iZ,iEq);
      iEq = iEq + 1;

      // east face
      calcSteadyStateEastCurrent(iR,iZ,iEq);
      iEq = iEq + 1;

      // north face
      if (iZ == 0)
      {
        calcSteadyStateNorthCurrent(iR,iZ,iEq);
        iEq = iEq + 1;
      } 

      // west face
      if (iR == 0)
      {
        // if on the boundary, assert boundary conditions
        calcSteadyStateWestCurrent(iR,iZ,iEq);
        iEq = iEq + 1;
      } 

    }
  }
};

//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate transient current on south face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcSouthCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->zNeutV(iZ+1,iR);
  double vPast = materials->oneGroupXS->zNeutVPast(iZ+1,iR);
  double sigT = materials->oneGroupXS->zSigTR(iZ+1,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,coeff,\
    zetaL,mgqdCurrent,mgqdNeutV; 
  double EzzC,EzzS,ErzW,ErzE;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zUp-zAvg;

  // get local Eddington factors
  EzzC = GGQD->Ezz(iZ,iR);
  EzzS = GGQD->EzzAxial(iZ+1,iR); 
  ErzW = GGQD->ErzRadial(iZ,iR);
  ErzE = GGQD->ErzRadial(iZ,iR+1);
  zetaL = materials->oneGroupXS->zZeta(iZ+1,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = 1/((1/(v*deltaT))+sigT); 

  C.coeffRef(iEq,indices[iSF]) -= coeff*EzzS/deltaZ;

  C.coeffRef(iEq,indices[iCF]) += coeff*EzzC/deltaZ;

  C.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzW/(rAvg*deltaR));

  C.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzE/(rAvg*deltaR);

  // Enforce zeta coefficient
  C.coeffRef(iEq,indices[iSF]) -= coeff*zetaL;

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      mgqdCurrent = GGQD->mgqd->QDSolve->currPast(indices[iSC]); 
      mgqdNeutV = materials->zNeutVel(iZ+1,iR,iGroup); 
      d(iEq) += (coeff/deltaT)*(mgqdCurrent/mgqdNeutV);
    }
  }
  else
  {
    d(iEq) = coeff*(currPast(indices[iSC])/(vPast*deltaT));
  }

};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate transient current on north face 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcNorthCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->zNeutV(iZ,iR);
  double vPast = materials->oneGroupXS->zNeutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->zSigTR(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,coeff,\
    zetaL,mgqdCurrent,mgqdNeutV;
  double EzzC,EzzN,ErzW,ErzE;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zAvg-zDown;

  // get local Eddington factors
  EzzC = GGQD->Ezz(iZ,iR);
  EzzN = GGQD->EzzAxial(iZ,iR); 
  ErzW = GGQD->ErzRadial(iZ,iR);
  ErzE = GGQD->ErzRadial(iZ,iR+1);

  zetaL = materials->oneGroupXS->zZeta(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = 1/((1/(v*deltaT))+sigT); 

  C.coeffRef(iEq,indices[iNF]) += coeff*EzzN/deltaZ;

  C.coeffRef(iEq,indices[iCF]) -= coeff*EzzC/deltaZ;

  C.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzW/(rAvg*deltaR));

  C.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzE/(rAvg*deltaR);

  // Enforce zeta coefficient
  C.coeffRef(iEq,indices[iNF]) -= coeff*zetaL;

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      mgqdCurrent = GGQD->mgqd->QDSolve->currPast(indices[iNC]); 
      mgqdNeutV = materials->zNeutVel(iZ,iR,iGroup); 
      d(iEq) += (coeff/deltaT)*(mgqdCurrent/mgqdNeutV);
    }
  }
  else
  {
    d(iEq) = coeff*(currPast(indices[iNC])/(vPast*deltaT));  
  }

};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate transient current on west face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcWestCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->rNeutV(iZ,iR);
  double vPast = materials->oneGroupXS->rNeutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->rSigTR(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,coeff;  
  double hCent,hDown,zetaL,mgqdCurrent,mgqdNeutV;
  double ErzN,ErzS,ErrC,ErrW;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rAvg-rDown; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg,iWF);
  hDown = calcIntegratingFactor(iR,iZ,rDown,iWF);

  // get local Eddington factors
  ErrC = GGQD->Err(iZ,iR);
  ErrW = GGQD->ErrRadial(iZ,iR); 
  ErzN = GGQD->ErzAxial(iZ,iR);
  ErzS = GGQD->ErzAxial(iZ+1,iR);
  zetaL = materials->oneGroupXS->rZeta(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = 1/((1/(v*deltaT))+sigT); 

  C.coeffRef(iEq,indices[iSF]) -= coeff*ErzS/deltaZ;

  C.coeffRef(iEq,indices[iNF]) += coeff*ErzN/deltaZ;

  C.coeffRef(iEq,indices[iCF]) -= coeff*hCent*ErrC/(hDown*deltaR);

  C.coeffRef(iEq,indices[iWF]) += coeff*hDown*ErrW/(hDown*deltaR);

  // Enforce zeta coefficient
  C.coeffRef(iEq,indices[iWF]) -= coeff*zetaL;

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      mgqdCurrent = GGQD->mgqd->QDSolve->currPast(indices[iWC]); 
      mgqdNeutV = materials->rNeutVel(iZ,iR,iGroup); 
      d(iEq) += (coeff/deltaT)*(mgqdCurrent/mgqdNeutV);
    }
  }
  else
  {
    d(iEq) = coeff*(currPast(indices[iWC])/(vPast*deltaT));
  }

};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate transient current on east face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcEastCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->rNeutV(iZ,iR+1);
  double vPast = materials->oneGroupXS->rNeutVPast(iZ,iR+1);
  double sigT = materials->oneGroupXS->rSigTR(iZ,iR+1);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,coeff;  
  double hCent,hUp,zetaL,mgqdCurrent,mgqdNeutV;
  double ErzN,ErzS,ErrC,ErrE;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rAvg; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg,iEF);
  hUp = calcIntegratingFactor(iR,iZ,rUp,iEF);

  // get local Eddington factors
  ErrC = GGQD->Err(iZ,iR);
  ErrE = GGQD->ErrRadial(iZ,iR+1); 
  ErzN = GGQD->ErzAxial(iZ,iR);
  ErzS = GGQD->ErzAxial(iZ+1,iR);

  zetaL = materials->oneGroupXS->rZeta(iZ,iR+1);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = 1/((1/(v*deltaT))+sigT); 

  C.coeffRef(iEq,indices[iSF]) -= coeff*ErzS/deltaZ;

  C.coeffRef(iEq,indices[iNF]) += coeff*ErzN/deltaZ;

  C.coeffRef(iEq,indices[iCF]) += coeff*hCent*ErrC/(hUp*deltaR);

  C.coeffRef(iEq,indices[iEF]) -= coeff*hUp*ErrE/(hUp*deltaR);

  // Enforce zeta coefficient
  C.coeffRef(iEq,indices[iEF]) -= coeff*zetaL;

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      mgqdCurrent = GGQD->mgqd->QDSolve->currPast(indices[iEC]); 
      mgqdNeutV = materials->rNeutVel(iZ,iR+1,iGroup); 
      d(iEq) += (coeff/deltaT)*(mgqdCurrent/mgqdNeutV);
    }
  }
  else
  {
    d(iEq) = coeff*(currPast(indices[iEC])/(vPast*deltaT));
  }

};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate steady state current on south face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcSteadyStateSouthCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->zNeutV(iZ+1,iR);
  double vPast = materials->oneGroupXS->zNeutVPast(iZ+1,iR);
  double sigT = materials->oneGroupXS->zSigTR(iZ+1,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,coeff,\
    zetaL,mgqdCurrent,mgqdNeutV; 
  double EzzC,EzzS,ErzW,ErzE;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zUp-zAvg;

  // get local Eddington factors
  EzzC = GGQD->Ezz(iZ,iR);
  EzzS = GGQD->EzzAxial(iZ+1,iR); 
  ErzW = GGQD->ErzRadial(iZ,iR);
  ErzE = GGQD->ErzRadial(iZ,iR+1);
  zetaL = materials->oneGroupXS->zZeta2(iZ+1,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = 1/(sigT); 

  C.coeffRef(iEq,indices[iSF]) -= coeff*EzzS/deltaZ;

  C.coeffRef(iEq,indices[iCF]) += coeff*EzzC/deltaZ;

  C.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzW/(rAvg*deltaR));

  C.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzE/(rAvg*deltaR);

  // Enforce zeta coefficient
  C.coeffRef(iEq,indices[iSF]) -= coeff*zetaL;

};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate steady state current on north face 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcSteadyStateNorthCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->zNeutV(iZ,iR);
  double vPast = materials->oneGroupXS->zNeutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->zSigTR(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,coeff,\
    zetaL,mgqdCurrent,mgqdNeutV;
  double EzzC,EzzN,ErzW,ErzE;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zAvg-zDown;

  // get local Eddington factors
  EzzC = GGQD->Ezz(iZ,iR);
  EzzN = GGQD->EzzAxial(iZ,iR); 
  ErzW = GGQD->ErzRadial(iZ,iR);
  ErzE = GGQD->ErzRadial(iZ,iR+1);

  zetaL = materials->oneGroupXS->zZeta2(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = 1/(sigT); 

  C.coeffRef(iEq,indices[iNF]) += coeff*EzzN/deltaZ;

  C.coeffRef(iEq,indices[iCF]) -= coeff*EzzC/deltaZ;

  C.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzW/(rAvg*deltaR));

  C.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzE/(rAvg*deltaR);

  // Enforce zeta coefficient
  C.coeffRef(iEq,indices[iNF]) -= coeff*zetaL;

};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate steady state current on west face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcSteadyStateWestCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->rNeutV(iZ,iR);
  double vPast = materials->oneGroupXS->rNeutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->rSigTR(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,coeff;  
  double hCent,hDown,zetaL,mgqdCurrent,mgqdNeutV;
  double ErzN,ErzS,ErrC,ErrW;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rAvg-rDown; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg,iWF);
  hDown = calcIntegratingFactor(iR,iZ,rDown,iWF);

  // get local Eddington factors
  ErrC = GGQD->Err(iZ,iR);
  ErrW = GGQD->ErrRadial(iZ,iR); 
  ErzN = GGQD->ErzAxial(iZ,iR);
  ErzS = GGQD->ErzAxial(iZ+1,iR);
  zetaL = materials->oneGroupXS->rZeta2(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = 1/(sigT); 

  C.coeffRef(iEq,indices[iSF]) -= coeff*ErzS/deltaZ;

  C.coeffRef(iEq,indices[iNF]) += coeff*ErzN/deltaZ;

  C.coeffRef(iEq,indices[iCF]) -= coeff*hCent*ErrC/(hDown*deltaR);

  C.coeffRef(iEq,indices[iWF]) += coeff*hDown*ErrW/(hDown*deltaR);

  // Enforce zeta coefficient
  C.coeffRef(iEq,indices[iWF]) -= coeff*zetaL;

};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate steady state current on east face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcSteadyStateEastCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->rNeutV(iZ,iR+1);
  double vPast = materials->oneGroupXS->rNeutVPast(iZ,iR+1);
  double sigT = materials->oneGroupXS->rSigTR(iZ,iR+1);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,coeff;  
  double hCent,hUp,zetaL,mgqdCurrent,mgqdNeutV;
  double ErzN,ErzS,ErrC,ErrE;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rAvg; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg,iEF);
  hUp = calcIntegratingFactor(iR,iZ,rUp,iEF);

  // get local Eddington factors
  ErrC = GGQD->Err(iZ,iR);
  ErrE = GGQD->ErrRadial(iZ,iR+1); 
  ErzN = GGQD->ErzAxial(iZ,iR);
  ErzS = GGQD->ErzAxial(iZ+1,iR);

  zetaL = materials->oneGroupXS->rZeta2(iZ,iR+1);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = 1/(sigT); 

  C.coeffRef(iEq,indices[iSF]) -= coeff*ErzS/deltaZ;

  C.coeffRef(iEq,indices[iNF]) += coeff*ErzN/deltaZ;

  C.coeffRef(iEq,indices[iCF]) += coeff*hCent*ErrC/(hUp*deltaR);

  C.coeffRef(iEq,indices[iEF]) -= coeff*hUp*ErrE/(hUp*deltaR);

  // Enforce zeta coefficient
  C.coeffRef(iEq,indices[iEF]) -= coeff*zetaL;

};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertNFluxBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);

  Atemp.insert(iEq,indices[iNF]) = 1.0;
  (*b)(iEq) = GGQD->nFluxBC(iR);
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSFluxBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);

  Atemp.insert(iEq,indices[iSF]) = 1.0;
  (*b)(iEq) = GGQD->sFluxBC(iR);
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the west face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertWFluxBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);

  Atemp.insert(iEq,indices[iWF]) = 1.0;
  (*b)(iEq) = GGQD->wFluxBC(iZ);
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertEFluxBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);

  Atemp.insert(iEq,indices[iEF]) = 1.0;
  (*b)(iEq) = GGQD->eFluxBC(iZ);
};
//==============================================================================

//==============================================================================
/// Assert the current boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertNCurrentBC(int iR,int iZ,int iEq)
{
  northCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert the current boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSCurrentBC(int iR,int iZ,int iEq)
{
  southCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert the current boundary condition on the west face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertWCurrentBC(int iR,int iZ,int iEq)
{
  westCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert the current boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertECurrentBC(int iR,int iZ,int iEq)
{
  eastCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert the steady state current boundary condition on the north face at 
/// location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateNCurrentBC(int iR,int iZ,int iEq)
{
  steadyStateNorthCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert the steady state current boundary condition on the south face at 
/// location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateSCurrentBC(int iR,int iZ,int iEq)
{
  steadyStateSouthCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert the steady state current boundary condition on the west face at 
/// location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateWCurrentBC(int iR,int iZ,int iEq)
{
  steadyStateWestCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert the steady state current boundary condition on the south face at 
/// location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateECurrentBC(int iR,int iZ,int iEq)
{
  steadyStateEastCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertNGoldinBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->nOutwardCurrToFluxRatioBC(iR);
  double inFluxWeightRatio = GGQD->nOutwardCurrToFluxRatioInwardWeightedBC(iR);
  double absCurrent = GGQD->nAbsCurrentBC(iR);
  double inwardCurrent = GGQD->nInwardCurrentBC(iR);
  double inwardFlux = GGQD->nInwardFluxBC(iR);

  northCurrent(1.0,iR,iZ,iEq);
  Atemp.coeffRef(iEq,indices[iNF]) -= ratio;
  (*b)(iEq) = (*b)(iEq) + (inwardCurrent-inFluxWeightRatio*inwardFlux);

  //northCurrent(1.0,iR,iZ,iEq,energyGroup,SGQD);
  //A.coeffRef(iEq,indices[iNF]) -= absCurrent/SGQD->nFluxBC(iR);
  //b(iEq) = b(iEq) + (2*inwardCurrent);

  //northCurrent(0.5,iR,iZ,iEq,energyGroup,SGQD);
  //southCurrent(0.5,iR,iZ,iEq,energyGroup,SGQD);
  //b(iEq) = b(iEq) + SGQD->nCurrentZBC(iR);

  //A.coeffRef(iEq,indices[iNF]) = 1.0;
  //b(iEq) = SGQD->nFluxBC(iR); 
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSGoldinBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->sOutwardCurrToFluxRatioBC(iR);
  double inFluxWeightRatio = GGQD->sOutwardCurrToFluxRatioInwardWeightedBC(iR);
  double absCurrent = GGQD->sAbsCurrentBC(iR);
  double inwardCurrent = GGQD->sInwardCurrentBC(iR);
  double inwardFlux = GGQD->sInwardFluxBC(iR);

  southCurrent(1.0,iR,iZ,iEq);
  Atemp.coeffRef(iEq,indices[iSF]) -= ratio;
  (*b)(iEq) = (*b)(iEq) + (inwardCurrent-inFluxWeightRatio*inwardFlux);

  //southCurrent(1.0,iR,iZ,iEq,energyGroup,SGQD);
  //A.coeffRef(iEq,indices[iSF]) -= absCurrent/SGQD->sFluxBC(iR);
  //b(iEq) = b(iEq) + (2*inwardCurrent);

  //southCurrent(0.5,iR,iZ,iEq,energyGroup,SGQD);
  //northCurrent(0.5,iR,iZ,iEq,energyGroup,SGQD);
  //b(iEq) = b(iEq) + SGQD->sCurrentZBC(iR);

  //A.coeffRef(iEq,indices[iSF]) = 1.0;
  //b(iEq) = SGQD->sFluxBC(iR); 
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertEGoldinBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->eOutwardCurrToFluxRatioBC(iZ);
  double inFluxWeightRatio = GGQD->eOutwardCurrToFluxRatioInwardWeightedBC(iZ);
  double absCurrent = GGQD->eAbsCurrentBC(iZ);
  double inwardCurrent = GGQD->eInwardCurrentBC(iZ);
  double inwardFlux = GGQD->eInwardFluxBC(iZ);

  eastCurrent(1.0,iR,iZ,iEq);
  Atemp.coeffRef(iEq,indices[iEF]) -= ratio;
  (*b)(iEq) = (*b)(iEq) + (inwardCurrent-inFluxWeightRatio*inwardFlux);

  //eastCurrent(1.0,iR,iZ,iEq,energyGroup,SGQD);
  //A.coeffRef(iEq,indices[iEF]) -= absCurrent/SGQD->eFluxBC(iZ);
  //b(iEq) = b(iEq) + (2*inwardCurrent);

  //eastCurrent(0.5,iR,iZ,iEq,energyGroup,SGQD);
  //westCurrent(0.5,iR,iZ,iEq,energyGroup,SGQD);
  //b(iEq) = b(iEq) + SGQD->eCurrentRBC(iZ);

  //A.coeffRef(iEq,indices[iEF]) = 1.0;
  //b(iEq) = SGQD->eFluxBC(iZ); 
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the north face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateNGoldinBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->nOutwardCurrToFluxRatioBC(iR);
  double inFluxWeightRatio = GGQD->nOutwardCurrToFluxRatioInwardWeightedBC(iR);
  double absCurrent = GGQD->nAbsCurrentBC(iR);
  double inwardCurrent = GGQD->nInwardCurrentBC(iR);
  double inwardFlux = GGQD->nInwardFluxBC(iR);

  steadyStateNorthCurrent(1.0,iR,iZ,iEq);
  Atemp.coeffRef(iEq,indices[iNF]) -= ratio;
  (*b)(iEq) = (*b)(iEq) + (inwardCurrent-inFluxWeightRatio*inwardFlux);

};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the south face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateSGoldinBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->sOutwardCurrToFluxRatioBC(iR);
  double inFluxWeightRatio = GGQD->sOutwardCurrToFluxRatioInwardWeightedBC(iR);
  double absCurrent = GGQD->sAbsCurrentBC(iR);
  double inwardCurrent = GGQD->sInwardCurrentBC(iR);
  double inwardFlux = GGQD->sInwardFluxBC(iR);

  steadyStateSouthCurrent(1.0,iR,iZ,iEq);
  Atemp.coeffRef(iEq,indices[iSF]) -= ratio;
  (*b)(iEq) = (*b)(iEq) + (inwardCurrent-inFluxWeightRatio*inwardFlux);

};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the east face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateEGoldinBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->eOutwardCurrToFluxRatioBC(iZ);
  double inFluxWeightRatio = GGQD->eOutwardCurrToFluxRatioInwardWeightedBC(iZ);
  double absCurrent = GGQD->eAbsCurrentBC(iZ);
  double inwardCurrent = GGQD->eInwardCurrentBC(iZ);
  double inwardFlux = GGQD->eInwardFluxBC(iZ);

  steadyStateEastCurrent(1.0,iR,iZ,iEq);
  Atemp.coeffRef(iEq,indices[iEF]) -= ratio;
  (*b)(iEq) = (*b)(iEq) + (inwardCurrent-inFluxWeightRatio*inwardFlux);

};
//==============================================================================

//==============================================================================
/// Assert boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertNBC(int iR,int iZ,int iEq)
{
  if (reflectingBCs)
    assertNCurrentBC(iR,iZ,iEq);
  else if (goldinBCs)
    assertNGoldinBC(iR,iZ,iEq);
  else
    assertNFluxBC(iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSBC(int iR,int iZ,int iEq)
{
  if (reflectingBCs)
    assertSCurrentBC(iR,iZ,iEq);
  else if (goldinBCs)
    assertSGoldinBC(iR,iZ,iEq);
  else
    assertSFluxBC(iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert boundary condition on the west face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertWBC(int iR,int iZ,int iEq)
{
  if (reflectingBCs or goldinBCs)
    assertWCurrentBC(iR,iZ,iEq);
  else
    // Can't think of a circumstance where there wouldn't be a reflecting BC at
    //   r = 0 
    //assertWFluxBC(iR,iZ,iEq);
    assertWCurrentBC(iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertEBC(int iR,int iZ,int iEq)
{
  if (reflectingBCs)
    assertECurrentBC(iR,iZ,iEq);
  else if (goldinBCs)
    assertEGoldinBC(iR,iZ,iEq);
  else
    assertEFluxBC(iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert steady state boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateNBC(int iR,int iZ,int iEq)
{
  if (reflectingBCs)
    assertSteadyStateNCurrentBC(iR,iZ,iEq);
  else if (goldinBCs)
    assertSteadyStateNGoldinBC(iR,iZ,iEq);
  else
    assertNFluxBC(iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert steady state boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateSBC(int iR,int iZ,int iEq)
{
  if (reflectingBCs)
    assertSteadyStateSCurrentBC(iR,iZ,iEq);
  else if (goldinBCs)
    assertSteadyStateSGoldinBC(iR,iZ,iEq);
  else
    assertSFluxBC(iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert steady state boundary condition on the west face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateWBC(int iR,int iZ,int iEq)
{
  if (reflectingBCs or goldinBCs)
    assertSteadyStateWCurrentBC(iR,iZ,iEq);
  else
    // Can't think of a circumstance where there wouldn't be a reflecting BC at
    //   r = 0 
    //assertWFluxBC(iR,iZ,iEq);
    assertSteadyStateWCurrentBC(iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert steady state boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertSteadyStateEBC(int iR,int iZ,int iEq)
{
  if (reflectingBCs)
    assertSteadyStateECurrentBC(iR,iZ,iEq);
  else if (goldinBCs)
    assertSteadyStateEGoldinBC(iR,iZ,iEq);
  else
    assertEFluxBC(iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Calculate multigroup source coefficient for cell at (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
double GreyGroupSolver::calcScatterAndFissionCoeff(int iR,int iZ)
{

  //double beta = GGQD->mpqd->mgdnp->beta;
  double localFissionSource,localSigF,localNu,localChiP,localSigS,\
    sourceCoefficient;

  localSigS = materials->oneGroupXS->sigS(iZ,iR);
  localFissionSource = materials->oneGroupXS->qdFluxCoeff(iZ,iR);
  sourceCoefficient = localSigS + localFissionSource;

  return sourceCoefficient;
};
//==============================================================================

////==============================================================================
///// Form a portion of the linear system that belongs to GGQD 
///// @param [in] iR radial index of cell
///// @param [in] iZ axial index of cell
///// @param [in] rEval radial location to evaluate integrating factor at
//double GreyGroupSolver::calcIntegratingFactor(int iR,int iZ,double rEval)	      
//{
//  double EzzL,ErrL,ErzL,G,rUp,rDown,rAvg,g0,g1,ratio,hEval;
//  int p;
//
//  // get local Eddington factors 
//  ErrL = GGQD->Err(iZ,iR);
//  EzzL = GGQD->Ezz(iZ,iR);
//  ErzL = GGQD->Erz(iZ,iR);
//
//  // evaluate G
//  //G = 1 + (ErrL+EzzL-1)/ErrL;
//  G = GGQD->G(iZ,iR);
//  G = 0;
//
//  // get boundaries of cell and calculate volume average 
//  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
//  rAvg = calcVolAvgR(rDown,rUp);
//
//  if (iR == 0){
//
//    // use a special expression for cells that share a boundary with
//    // the z-axis
//    p = 2;
//    ratio = (pow(rUp,p+1)-pow(rAvg,p+1))/(pow(rAvg,p)-pow(rUp,p));
//    g1 = G/(pow(rAvg,p) * (rAvg + ratio));
//    g0 = g1 * ratio;
//
//    hEval = exp((g0*pow(rEval,p)/p)+g1*(pow(rEval,p+1))/(p+1));
//
//  } else {
//
//    // use the typical expression
//    hEval = pow(rEval,G);
//
//  }
//
//  return hEval;
//
//};
//
////==============================================================================

//==============================================================================
/// Form a portion of the linear system that belongs to GGQD 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] rEval radial location to evaluate integrating factor at
double GreyGroupSolver::calcIntegratingFactor(int iR,int iZ,double rEval,int iLoc) 
{
  double g0,g1,hEval,G;
  int p = 2;

  if (iR == 0)
  {
    // use a special expression for cells that share a boundary with
    // the z-axis
    g1 = GGQD->g1(iZ);
    g0 = GGQD->g0(iZ);

    hEval = exp((g0*pow(rEval,p)/p)+g1*(pow(rEval,p+1))/(p+1));
  } 
  else {

    // use the typical expressions depending on location to get G for
    if (iLoc == iWF)
    {
      G = GGQD->GL(iZ,iR);
      hEval = pow(rEval,G);
    }
    else if (iLoc == iEF)
    {
      G = GGQD->GR(iZ,iR);
      hEval = pow(rEval,G);
    }
  }

  return hEval;

};

//==============================================================================


//==============================================================================
/// Return global index of south face current at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for south face current in cell at (iR,iZ)
vector<int> GreyGroupSolver::getIndices(int iR,int iZ)
{

  vector<int> indices,oneGroupIndices;

  // Get indices for a single energy group 
  oneGroupIndices = mesh->getQDCellIndices(iR,iZ);

  // Offset by specified energy group
  indices.push_back(oneGroupIndices[iCF] + GGQD->indexOffset);
  indices.push_back(oneGroupIndices[iWF] + GGQD->indexOffset);
  indices.push_back(oneGroupIndices[iEF] + GGQD->indexOffset);
  indices.push_back(oneGroupIndices[iNF] + GGQD->indexOffset);
  indices.push_back(oneGroupIndices[iSF] + GGQD->indexOffset);
  indices.push_back(oneGroupIndices[iWC]);
  indices.push_back(oneGroupIndices[iEC]);
  indices.push_back(oneGroupIndices[iNC]);
  indices.push_back(oneGroupIndices[iSC]);

  return indices;
};

//==============================================================================


//==============================================================================
/// Return geometry parameters for the cell located at (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] gParams vector containing volume and surfaces areas of the 
///   west, east, north, and south faces, in that order.
vector<double> GreyGroupSolver::calcGeoParams(int iR,int iZ)
{
  vector<double> gParams;
  double rDown,rUp,zDown,zUp,volume,wFaceSA,eFaceSA,nFaceSA,sFaceSA;

  // get boundaries of this cell
  rDown = mesh->rCornerEdge(iR); rUp = mesh->rCornerEdge(iR+1);
  zDown = mesh->zCornerEdge(iR); zUp = mesh->zCornerEdge(iR+1);

  // calculate geometry parameters
  volume = M_PI*(rUp*rUp-rDown*rDown)*(zUp-zDown);
  nFaceSA = M_PI*(rUp*rUp-rDown*rDown); 
  sFaceSA = nFaceSA;
  eFaceSA = 2*M_PI*rUp*(zUp-zDown);
  wFaceSA = 2*M_PI*rDown*(zUp-zDown);

  // add parameters to the vector and return it
  gParams.push_back(volume);
  gParams.push_back(wFaceSA);
  gParams.push_back(eFaceSA);
  gParams.push_back(nFaceSA);
  gParams.push_back(sFaceSA);

  return gParams;
};
//==============================================================================

//==============================================================================
/// Return volume-averaged radial coordinate for cell with boundaries rDown and
///   rUp
/// @param [in] rDown location of left cell edge
/// @param [in] rUp location of right cell edge
/// @param [out] volAvgR volume-averaged radial coordinate 
double GreyGroupSolver::calcVolAvgR(double rDown,double rUp)
{
  // calculate volume-averaged radius
  double volAvgR = (2.0/3.0)*(pow(rUp,3) - pow(rDown,3))/(pow(rUp,2)-pow(rDown,2));

  return volAvgR;
};
//==============================================================================

//==============================================================================
/// Get Err on west interface 
///
/// @param [in] iZ axial coordinate 
/// @param [in] iR radial coordinate 
/// @param [out] Eddington value to use at interface 
double GreyGroupSolver::getWestErr(int iZ,int iR)
{

  double volLeft,volRight,totalVol;

  // Return the harmonic average of the Eddingtons on either side of the west
  // interface, unless on a boundary. In that case, simply return the cell
  // center Eddington.
  if (iR == 0)
    return GGQD->Err(iZ,iR);
  else
  { 
    volLeft = calcGeoParams(iR-1,iZ)[iCF];
    volRight = calcGeoParams(iR,iZ)[iCF];
    totalVol = volLeft + volRight;
    return (totalVol)/(volLeft/GGQD->Err(iZ,iR-1) + volRight/GGQD->Err(iZ,iR));
  }

};
//==============================================================================

//==============================================================================
/// Get Erz on west interface 
///
/// @param [in] iZ axial coordinate 
/// @param [in] iR radial coordinate 
/// @param [out] Eddington value to use at interface 
double GreyGroupSolver::getWestErz(int iZ,int iR)
{

  double volLeft,volRight,totalVol;

  // Return the harmonic average of the Eddingtons on either side of the west
  // interface, unless on a boundary. In that case, simply return the cell
  // center Eddington.
  if (iR == 0)
    return GGQD->Erz(iZ,iR);
  else
  { 
    volLeft = calcGeoParams(iR-1,iZ)[iCF];
    volRight = calcGeoParams(iR,iZ)[iCF];
    totalVol = volLeft + volRight;
    return (totalVol)/(volLeft/GGQD->Erz(iZ,iR-1) + volRight/GGQD->Erz(iZ,iR));
  }

};
//==============================================================================

//==============================================================================
/// Get Err on east interface 
///
/// @param [in] iZ axial coordinate 
/// @param [in] iR radial coordinate 
/// @param [out] Eddington value to use at interface 
double GreyGroupSolver::getEastErr(int iZ,int iR)
{

  double volLeft,volRight,totalVol;

  // Return the harmonic average of the Eddingtons on either side of the east
  // interface, unless on a boundary. In that case, simply return the cell
  // center Eddington.
  if (iR == mesh->nR - 1)
    return GGQD->Err(iZ,iR);
  else
  { 
    volLeft = calcGeoParams(iR,iZ)[iCF];
    volRight = calcGeoParams(iR+1,iZ)[iCF];
    totalVol = volLeft + volRight;
    return (totalVol)/(volLeft/GGQD->Err(iZ,iR) + volRight/GGQD->Err(iZ,iR+1));
  }

};
//==============================================================================

//==============================================================================
/// Get Erz on east interface 
///
/// @param [in] iZ axial coordinate 
/// @param [in] iR radial coordinate 
/// @param [out] Eddington value to use at interface 
double GreyGroupSolver::getEastErz(int iZ,int iR)
{

  double volLeft,volRight,totalVol;

  // Return the harmonic average of the Eddingtons on either side of the east
  // interface, unless on a boundary. In that case, simply return the cell
  // center Eddington.
  if (iR == mesh->nR - 1)
    return GGQD->Erz(iZ,iR);
  else
  { 
    volLeft = calcGeoParams(iR,iZ)[iCF];
    volRight = calcGeoParams(iR+1,iZ)[iCF];
    totalVol = volLeft + volRight;
    return (totalVol)/(volLeft/GGQD->Erz(iZ,iR) + volRight/GGQD->Erz(iZ,iR+1));
  }

};
//==============================================================================

//==============================================================================
/// Get Ezz on north interface 
///
/// @param [in] iZ axial coordinate 
/// @param [in] iR radial coordinate 
/// @param [out] Eddington value to use at interface 
double GreyGroupSolver::getNorthEzz(int iZ,int iR)
{

  double volDown,volUp,totalVol;

  // Return the harmonic average of the Eddingtons on either side of the north
  // interface, unless on a boundary. In that case, simply return the cell
  // center Eddington.
  if (iZ == 0)
    return GGQD->Ezz(iZ,iR);
  else
  { 
    volDown = calcGeoParams(iR,iZ-1)[iCF];
    volUp = calcGeoParams(iR,iZ)[iCF];
    totalVol = volUp + volDown;
    return (totalVol)/(volDown/GGQD->Ezz(iZ-1,iR) + volUp/GGQD->Ezz(iZ,iR));
  }

};
//==============================================================================

//==============================================================================
/// Get Erz on north interface 
///
/// @param [in] iZ axial coordinate 
/// @param [in] iR radial coordinate 
/// @param [out] Eddington value to use at interface 
double GreyGroupSolver::getNorthErz(int iZ,int iR)
{

  double volDown,volUp,totalVol;

  // Return the harmonic average of the Eddingtons on either side of the north
  // interface, unless on a boundary. In that case, simply return the cell
  // center Eddington.
  if (iZ == 0)
    return GGQD->Erz(iZ,iR);
  else
  { 
    volDown = calcGeoParams(iR,iZ-1)[iCF];
    volUp = calcGeoParams(iR,iZ)[iCF];
    totalVol = volUp + volDown;
    return (totalVol)/(volDown/GGQD->Erz(iZ-1,iR) + volUp/GGQD->Erz(iZ,iR));
  }

};
//==============================================================================

//==============================================================================
/// Get Ezz on south interface 
///
/// @param [in] iZ axial coordinate 
/// @param [in] iR radial coordinate 
/// @param [out] Eddington value to use at interface 
double GreyGroupSolver::getSouthEzz(int iZ,int iR)
{

  double volDown,volUp,totalVol;

  // Return the harmonic average of the Eddingtons on either side of the south
  // interface, unless on a boundary. In that case, simply return the cell
  // center Eddington.
  if (iZ == mesh->nZ-1)
    return GGQD->Ezz(iZ,iR);
  else
  { 
    volDown = calcGeoParams(iR,iZ)[iCF];
    volUp = calcGeoParams(iR,iZ+1)[iCF];
    totalVol = volUp + volDown;
    return (totalVol)/(volDown/GGQD->Ezz(iZ,iR) + volUp/GGQD->Ezz(iZ+1,iR));
  }

};
//==============================================================================

//==============================================================================
/// Get Erz on south interface 
///
/// @param [in] iZ axial coordinate 
/// @param [in] iR radial coordinate 
/// @param [out] Eddington value to use at interface 
double GreyGroupSolver::getSouthErz(int iZ,int iR)
{

  double volDown,volUp,totalVol;

  // Return the harmonic average of the Eddingtons on either side of the south
  // interface, unless on a boundary. In that case, simply return the cell
  // center Eddington.
  if (iZ == mesh->nZ-1)
    return GGQD->Erz(iZ,iR);
  else
  { 
    volDown = calcGeoParams(iR,iZ)[iCF];
    volUp = calcGeoParams(iR,iZ+1)[iCF];
    totalVol = volUp + volDown;
    return (totalVol)/(volDown/GGQD->Erz(iZ,iR) + volUp/GGQD->Erz(iZ+1,iR));
  }

};
//==============================================================================

//==============================================================================
/// Extract cell average values from solution vector and store
void GreyGroupSolver::getFlux()
{
  vector<int> indices;

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      indices = getIndices(iR,iZ);

      // Read fluxes into recirculation back calc vector
      xFlux(indices[iCF]) = (*x)(indices[iCF]);
      xFlux(indices[iWF]) = (*x)(indices[iWF]);
      xFlux(indices[iEF]) = (*x)(indices[iEF]);
      xFlux(indices[iNF]) = (*x)(indices[iNF]);
      xFlux(indices[iSF]) = (*x)(indices[iSF]);

      // Read fluxes into GGQD object
      GGQD->sFlux(iZ,iR) = xFlux(indices[iCF]);
      GGQD->sFluxR(iZ,iR) = xFlux(indices[iWF]);
      GGQD->sFluxR(iZ,iR+1) = xFlux(indices[iEF]);
      GGQD->sFluxZ(iZ,iR) = xFlux(indices[iNF]);
      GGQD->sFluxZ(iZ+1,iR) = xFlux(indices[iSF]);

    }
  } 

};
//==============================================================================

//==============================================================================
/// Extract cell average values from solution vector and store
void GreyGroupSolver::getCurrent()
{
  vector<int> indices;

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      indices = getIndices(iR,iZ);

      // Get cell currents 
      GGQD->currentR(iZ,iR) = currPast(indices[iWC]);
      GGQD->currentR(iZ,iR+1) = currPast(indices[iEC]);
      GGQD->currentZ(iZ,iR) = currPast(indices[iNC]);
      GGQD->currentZ(iZ+1,iR) = currPast(indices[iSC]);

    }
  } 

};
//==============================================================================

//==============================================================================
/// Map values from 2D matrix to 1D solution vector
void GreyGroupSolver::setFlux()
{
  vector<int> indices;

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      indices = getIndices(iR,iZ);
      (*xPast)(indices[iCF]) = GGQD->sFlux(iZ,iR);

      (*xPast)(indices[iWF]) = GGQD->sFluxR(iZ,iR);
      (*xPast)(indices[iEF]) = GGQD->sFluxR(iZ,iR+1);
      (*xPast)(indices[iNF]) = GGQD->sFluxZ(iZ,iR);
      (*xPast)(indices[iSF]) = GGQD->sFluxZ(iZ+1,iR);

    }
  }  

};
//==============================================================================

//==============================================================================
/// Extract flux values from GGQD object, store to solution vector, and return 
/// @param [out] solVector flux values mapped to the 1D vector
Eigen::VectorXd GreyGroupSolver::getFluxSolutionVector()
{
  Eigen::VectorXd solVector(nUnknowns);
  solVector.setZero();
  vector<int> indices;

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      indices = getIndices(iR,iZ);
      solVector(indices[iCF]) = GGQD->sFlux(iZ,iR);

      solVector(indices[iWF]) = GGQD->sFluxR(iZ,iR);
      solVector(indices[iEF]) = GGQD->sFluxR(iZ,iR+1);
      solVector(indices[iNF]) = GGQD->sFluxZ(iZ,iR);
      solVector(indices[iSF]) = GGQD->sFluxZ(iZ+1,iR);

    }
  }  

  return solVector;
};
//==============================================================================

//==============================================================================
/// Extract current values from GGQD object, store to solution vector, and 
/// return 
/// @param [out] solVector current values mapped to the 1D vector
Eigen::VectorXd GreyGroupSolver::getCurrentSolutionVector()
{
  Eigen::VectorXd solVector(nCurrentUnknowns);
  solVector.setZero();
  vector<int> indices;

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      indices = getIndices(iR,iZ);

      solVector(indices[iWC]) = GGQD->currentR(iZ,iR);
      solVector(indices[iEC]) = GGQD->currentR(iZ,iR+1);
      solVector(indices[iNC]) = GGQD->currentZ(iZ,iR);
      solVector(indices[iSC]) = GGQD->currentZ(iZ+1,iR);
    }
  }  

  return solVector;
};
//==============================================================================

//==============================================================================
/// Assign pointers to linear system components 
///
void GreyGroupSolver::assignPointers(Eigen::SparseMatrix<double,Eigen::RowMajor> * myA,\
    Eigen::VectorXd * myx,\
    Eigen::VectorXd * myxpast,\
    Eigen::VectorXd * myb)
{

  A = myA;
  x = myx;
  xPast = myxpast;
  b = myb;

};
//==============================================================================


//==============================================================================
/// Check for optional inputs of relevance to this object
void GreyGroupSolver::checkOptionalParams()
{
  string boundaryType;

  // check for optional parameters specified in input file

  if ((*input)["parameters"]["solve type"])
  {

    boundaryType=(*input)["parameters"]["solve type"].as<string>();
    if (boundaryType == "TQD") goldinBCs = true;

  }

  if ((*input)["parameters"]["mgqd-bcs"])
  {

    boundaryType=(*input)["parameters"]["mgqd-bcs"].as<string>();

    if (boundaryType == "reflective" or boundaryType == "REFLECTIVE"\
        or boundaryType == "Reflective")
      reflectingBCs = true;
    else if (boundaryType == "goldin" or boundaryType == "GOLDIN" \
        or boundaryType == "Goldin")
      goldinBCs = true;

  }
}
//==============================================================================
