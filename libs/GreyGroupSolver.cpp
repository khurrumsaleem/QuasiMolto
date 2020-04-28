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

  checkOptionalParams();
};

//==============================================================================

//==============================================================================
/// Form a portion of the linear system that belongs to SGQD 
/// @param [in] GGQD quasidiffusion energy group to build portion of linear 
///   for
void GreyGroupSolver::formLinearSystem(int iEq)	      
{

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {

      // apply zeroth moment equation
      assertZerothMoment(iR,iZ,iEq);
      iEq = iEq + 1;

      // north face
      if (iZ == 0)
      {
        // if on the boundary, assert boundary conditions
        assertNBC(iR,iZ,iEq);
        iEq = iEq + 1;
      } 

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

      // west face
      if (iR == 0)
      {
        // if on the boundary, assert boundary conditions
        assertWBC(iR,iZ,iEq);
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
    }
  }
};

//==============================================================================

//==============================================================================
/// Compute currents from flux values in x
void GreyGroupSolver::backCalculateCurrent()
{
  currPast = d + C*(xFlux); 
}
//==============================================================================


//==============================================================================
/// Assert the zeroth moment equation for cell (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertZerothMoment(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double groupSourceCoeff;

  indices = getIndices(iR,iZ);
  
  // Scatter and fission source term
  groupSourceCoeff = calcScatterAndFissionCoeff(iZ,iR);
  A->insert(iEq,indices[iCF]) = -geoParams[iCF] * groupSourceCoeff;

  // DNP source term
  GGQD->mpqd->dnpSource(iZ,iR,iEq,-geoParams[iCF]);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  A->coeffRef(iEq,indices[iCF]) += geoParams[iCF] * ((1/(v*deltaT)) + sigT);

  westCurrent(-geoParams[iWF],iR,iZ,iEq);
  
  eastCurrent(geoParams[iEF],iR,iZ,iEq);

  northCurrent(-geoParams[iNF],iR,iZ,iEq);

  southCurrent(geoParams[iSF],iR,iZ,iEq);

  // formulate RHS entry
  (*b)(iEq) = (*b)(iEq) + geoParams[iCF]*\
    ( ((*xPast)(indices[iCF])/(v*deltaT)) + GGQD->q(iZ,iR));
};
//==============================================================================

//==============================================================================
/// Apply radial boundary for cell (iR,iZ)
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
/// Apply axial boundary for cell (iR,iZ)
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
/// Enforce coefficients for current on south face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::southCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zUp-zAvg;

  // get local Eddington factors
  ErrL = GGQD->Err(iZ,iR);
  EzzL = GGQD->Ezz(iZ,iR);
  ErzL = GGQD->Erz(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = coeff/((1/(v*deltaT))+sigT); 

  A->coeffRef(iEq,indices[iSF]) -= coeff*EzzL/deltaZ;

  A->coeffRef(iEq,indices[iCF]) += coeff*EzzL/deltaZ;

  A->coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzL/(rAvg*deltaR));

  A->coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzL/(rAvg*deltaR);
  
  // formulate RHS entry
  (*b)(iEq) = (*b)(iEq) - coeff*(currPast(indices[iSC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for current on north face 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::northCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zAvg-zDown;

  // get local Eddington factors
  ErrL = GGQD->Err(iZ,iR);
  EzzL = GGQD->Ezz(iZ,iR);
  ErzL = GGQD->Erz(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);
  
  coeff = coeff/((1/(v*deltaT))+sigT); 

  A->coeffRef(iEq,indices[iNF]) += coeff*EzzL/deltaZ;

  A->coeffRef(iEq,indices[iCF]) -= coeff*EzzL/deltaZ;

  A->coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzL/(rAvg*deltaR));

  A->coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzL/(rAvg*deltaR);

  // formulate RHS entry
  (*b)(iEq) = (*b)(iEq) - coeff*(currPast(indices[iNC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for current on west face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::westCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL;  
  double hCent,hDown;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rAvg-rDown; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg);
  hDown = calcIntegratingFactor(iR,iZ,rDown);

  // get local Eddington factors
  ErrL = GGQD->Err(iZ,iR);
  EzzL = GGQD->Ezz(iZ,iR);
  ErzL = GGQD->Erz(iZ,iR);
  
  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);
  
  coeff = coeff/((1/(v*deltaT))+sigT); 

  A->coeffRef(iEq,indices[iSF]) -= coeff*ErzL/deltaZ;

  A->coeffRef(iEq,indices[iNF]) += coeff*ErzL/deltaZ;

  A->coeffRef(iEq,indices[iCF]) -= coeff*hCent*ErrL/(hDown*deltaR);

  A->coeffRef(iEq,indices[iWF]) += coeff*hDown*ErrL/(hDown*deltaR);

  // formulate RHS entry
  (*b)(iEq) = (*b)(iEq) - coeff*(currPast(indices[iWC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for current on east face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::eastCurrent(double coeff,int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL;  
  double hCent,hUp;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rAvg; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg);
  hUp = calcIntegratingFactor(iR,iZ,rUp);

  // get local Eddington factors
  ErrL = GGQD->Err(iZ,iR);
  EzzL = GGQD->Ezz(iZ,iR);
  ErzL = GGQD->Erz(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);
  
  coeff = coeff/((1/(v*deltaT))+sigT); 

  A->coeffRef(iEq,indices[iSF]) -= coeff*ErzL/deltaZ;

  A->coeffRef(iEq,indices[iNF]) += coeff*ErzL/deltaZ;

  A->coeffRef(iEq,indices[iCF]) += coeff*hCent*ErrL/(hUp*deltaR);

  A->coeffRef(iEq,indices[iEF]) -= coeff*hUp*ErrL/(hUp*deltaR);
  
  // formulate RHS entry
  (*b)(iEq) = (*b)(iEq) - coeff*(currPast(indices[iEC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Form a portion of the current back calc linear system that belongs to GGQD 
/// @param [in] GGQD quasidiffusion energy group to build portion of linear 
///   for
void GreyGroupSolver::formBackCalcSystem()	      
{
  int iEq = GGQD->indexOffset;

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
/// Enforce coefficients to calculate current on south face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcSouthCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL,coeff;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zUp-zAvg;

  // get local Eddington factors
  ErrL = GGQD->Err(iZ,iR);
  EzzL = GGQD->Ezz(iZ,iR);
  ErzL = GGQD->Erz(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  coeff = 1/((1/(v*deltaT))+sigT); 

  C.coeffRef(iEq,indices[iSF]) -= coeff*EzzL/deltaZ;

  C.coeffRef(iEq,indices[iCF]) += coeff*EzzL/deltaZ;

  C.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzL/(rAvg*deltaR));

  C.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzL/(rAvg*deltaR);
  
  // formulate RHS entry
  d(iEq) = coeff*(currPast(indices[iSC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate current on north face 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcNorthCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL,coeff;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zAvg-zDown;

  // get local Eddington factors
  ErrL = GGQD->Err(iZ,iR);
  EzzL = GGQD->Ezz(iZ,iR);
  ErzL = GGQD->Erz(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);
  
  coeff = 1/((1/(v*deltaT))+sigT); 

  C.coeffRef(iEq,indices[iNF]) += coeff*EzzL/deltaZ;

  C.coeffRef(iEq,indices[iCF]) -= coeff*EzzL/deltaZ;

  C.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzL/(rAvg*deltaR));

  C.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzL/(rAvg*deltaR);

  // formulate RHS entry
  d(iEq) = coeff*(currPast(indices[iNC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate current on west face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcWestCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL,coeff;  
  double hCent,hDown;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rAvg-rDown; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg);
  hDown = calcIntegratingFactor(iR,iZ,rDown);

  // get local Eddington factors
  ErrL = GGQD->Err(iZ,iR);
  EzzL = GGQD->Ezz(iZ,iR);
  ErzL = GGQD->Erz(iZ,iR);
  
  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);
  
  coeff = 1/((1/(v*deltaT))+sigT); 

  C.coeffRef(iEq,indices[iSF]) -= coeff*ErzL/deltaZ;

  C.coeffRef(iEq,indices[iNF]) += coeff*ErzL/deltaZ;

  C.coeffRef(iEq,indices[iCF]) -= coeff*hCent*ErrL/(hDown*deltaR);

  C.coeffRef(iEq,indices[iWF]) += coeff*hDown*ErrL/(hDown*deltaR);

  // formulate RHS entry
  d(iEq) = coeff*(currPast(indices[iWC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate current on east face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::calcEastCurrent(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL,coeff;  
  double hCent,hUp;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rAvg; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg);
  hUp = calcIntegratingFactor(iR,iZ,rUp);

  // get local Eddington factors
  ErrL = GGQD->Err(iZ,iR);
  EzzL = GGQD->Ezz(iZ,iR);
  ErzL = GGQD->Erz(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);
  
  coeff = 1/((1/(v*deltaT))+sigT); 

  C.coeffRef(iEq,indices[iSF]) -= coeff*ErzL/deltaZ;

  C.coeffRef(iEq,indices[iNF]) += coeff*ErzL/deltaZ;

  C.coeffRef(iEq,indices[iCF]) += coeff*hCent*ErrL/(hUp*deltaR);

  C.coeffRef(iEq,indices[iEF]) -= coeff*hUp*ErrL/(hUp*deltaR);
  
  // formulate RHS entry
  d(iEq) = coeff*(currPast(indices[iEC])/(v*deltaT));
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
  
  A->insert(iEq,indices[iNF]) = 1.0;
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
  
  A->insert(iEq,indices[iSF]) = 1.0;
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
  
  A->insert(iEq,indices[iWF]) = 1.0;
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
  
  A->insert(iEq,indices[iEF]) = 1.0;
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
/// Assert Gol'din's boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertNGoldinBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->nOutwardCurrToFluxRatioBC(iR);
  double absCurrent = GGQD->nAbsCurrentBC(iR);
  double inwardCurrent = GGQD->nInwardCurrentBC(iR);
  double inwardFlux = GGQD->nInwardFluxBC(iR);

  northCurrent(1.0,iR,iZ,iEq);
  A->coeffRef(iEq,indices[iNF]) -= ratio;
  (*b)(iEq) = (*b)(iEq) + (inwardCurrent-ratio*inwardFlux);
  
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
  double absCurrent = GGQD->sAbsCurrentBC(iR);
  double inwardCurrent = GGQD->sInwardCurrentBC(iR);
  double inwardFlux = GGQD->sInwardFluxBC(iR);

  southCurrent(1.0,iR,iZ,iEq);
  A->coeffRef(iEq,indices[iSF]) -= ratio;
  (*b)(iEq) = (*b)(iEq) + (inwardCurrent-ratio*inwardFlux);

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
  double absCurrent = GGQD->eAbsCurrentBC(iZ);
  double inwardCurrent = GGQD->eInwardCurrentBC(iZ);
  double inwardFlux = GGQD->eInwardFluxBC(iZ);

  eastCurrent(1.0,iR,iZ,iEq);
  A->coeffRef(iEq,indices[iEF]) -= ratio;
  (*b)(iEq) = (*b)(iEq) + (inwardCurrent-ratio*inwardFlux);
  
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
    assertWFluxBC(iR,iZ,iEq);
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
/// Calculate multigroup source coefficient for cell at (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
double GreyGroupSolver::calcScatterAndFissionCoeff(int iR,int iZ)
{

  double beta = GGQD->mpqd->mgdnp->beta;
  double localSigF,localNu,localChiP,localSigS,sourceCoefficient;

  localSigF = materials->oneGroupXS->sigF(iZ,iR);
  localNu = materials->oneGroupXS->nu(iZ,iR);
  localSigS = materials->oneGroupXS->sigS(iZ,iR);
  sourceCoefficient = localSigS + (1-beta)*localNu*localSigF;

  return sourceCoefficient;
};
//==============================================================================

//==============================================================================
/// Form a portion of the linear system that belongs to GGQD 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] rEval radial location to evaluate integrating factor at
double GreyGroupSolver::calcIntegratingFactor(int iR,int iZ,double rEval)	      
{
  double EzzL,ErrL,ErzL,G,rUp,rDown,rAvg,g0,g1,ratio,hEval;
  int p;
  
  // get local Eddington factors 
  ErrL = GGQD->Err(iZ,iR);
  EzzL = GGQD->Ezz(iZ,iR);
  ErzL = GGQD->Erz(iZ,iR);

  // evaluate G
  G = 1 + (ErrL+EzzL-1)/ErrL;

  // get boundaries of cell and calculate volume average 
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  rAvg = calcVolAvgR(rDown,rUp);

  if (iR == 0){

    // use a special expression for cells that share a boundary with
    // the z-axis
    p = 2;
    ratio = (pow(rUp,p+1)-pow(rAvg,p+1))/(pow(rAvg,p)-pow(rUp,p));
    g1 = G/(pow(rAvg,p) * (rAvg + ratio));
    g0 = g1 * ratio;

    hEval = exp((g0*pow(rEval,p)/p)+g1*(pow(rEval,p+1))/(p+1));

  } else {

    // use the typical expression
    hEval = pow(rEval,G);

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
      xFlux(indices[iCF]) = (*x)(indices[iCF]);
      GGQD->sFlux(iZ,iR) = xFlux(indices[iCF]);

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
void GreyGroupSolver::assignPointers(Eigen::SparseMatrix<double> * myA,\
  Eigen::VectorXd * myx,\
  Eigen::VectorXd * myb)
{

  A = myA;
  x = myx;
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
