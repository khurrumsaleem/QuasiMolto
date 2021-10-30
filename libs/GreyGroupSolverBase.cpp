// File: GreyGroupSolverBase.cpp
// Purpose: Solve RZ quasidiffusion equations  
// Date: February 05, 2020

#include "GreyGroupSolverBase.h"
#include "GreyGroupQD.h"

using namespace std; 

//==============================================================================
/// GreyGroupSolver object constructor
///
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myInput YAML input object for the simulation
/// @param [in] myMPQD multiphysics coupled QD object for the simulation
GreyGroupSolverBase::GreyGroupSolverBase(GreyGroupQD * myGGQD,\
    Mesh * myMesh,\
    Materials * myMaterials,\
    YAML::Node * myInput)	      
{
  // PETSc object to broadcast variables 
  VecScatter     ctx;

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

  /* Initialize PETSc variables */
  // Multiphysics system variables 
  initPETScRectMat(&C,nCurrentUnknowns,nUnknowns,20);
  initPETScVec(&currPast,nCurrentUnknowns);
  initPETScVec(&d,nCurrentUnknowns);
  initPETScVec(&xFlux,nUnknowns);

  // Broadcast currPast
  VecScatterCreateToAll(currPast,&ctx,&(currPastSeq));
  VecScatterBegin(ctx,currPast,currPastSeq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,currPast,currPastSeq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);

  // Set number of unknowns if GreyGroupQD object
  GGQD->nUnknowns = nUnknowns;
  GGQD->nCurrentUnknowns = nCurrentUnknowns;

  checkOptionalParams();
};

//==============================================================================

//==============================================================================
/// GreyGroupSolver copy constructor
///
/// @param [in] myMPQD multiphysics coupled QD object for the simulation
GreyGroupSolverBase::GreyGroupSolverBase(const GreyGroupSolverBase &solver)
  : GGQD(solver.GGQD) 
  , mesh(solver.mesh) 
  , input(solver.input) 
  , materials(solver.materials) 
  , nR(solver.nR) 
  , nZ(solver.nZ) 
  , nUnknowns(solver.nUnknowns) 
  , nCurrentUnknowns(solver.nCurrentUnknowns) 
  , C(solver.C)
  , currPast(solver.currPast)
  , d(solver.d)
  , xFlux(solver.xFlux)
  , currPastSeq(solver.currPastSeq)
{
  
  checkOptionalParams();

};

//==============================================================================


////==============================================================================
/// Form a portion of the linear system  

void GreyGroupSolverBase::formLinearSystem()	      
{

  int iEq = GGQD->indexOffset;

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
};

//==============================================================================

////==============================================================================
/// Apply radial boundary for cell (iR,iZ)
///
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolverBase::applyRadialBoundary(int iR,int iZ,int iEq)
{
  eastCurrent(1,iR,iZ,iEq);
  westCurrent(-1,iR+1,iZ,iEq);
}
//==============================================================================

//==============================================================================
/// Apply axial boundary for cell (iR,iZ)
///
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolverBase::applyAxialBoundary(int iR,int iZ,int iEq)
{
  northCurrent(1,iR,iZ+1,iEq);
  southCurrent(-1,iR,iZ,iEq);
}
//==============================================================================

////==============================================================================
/// Form a portion of the current back calc linear system that belongs to GGQD 
///
int GreyGroupSolverBase::formBackCalcSystem()	      
{
  int iEq = GGQD->indexOffset;
  PetscErrorCode ierr;

  // Reset linear system
  MatZeroEntries(C);
  VecZeroEntries(d);

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

  /* Finalize assembly for C and d */
  ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  
  ierr = VecAssemblyBegin(d);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(d);CHKERRQ(ierr);

  return ierr;

};

//==============================================================================

////==============================================================================
/// Compute currents from flux values in x
///
int GreyGroupSolverBase::backCalculateCurrent()
{
  // PETSc object to broadcast variables 
  VecScatter     ctx;
  PetscErrorCode ierr;

  ierr = MatMultAdd(C,xFlux,d,currPast);

  // Broadcast currPast
  VecDestroy(&(currPastSeq));
  VecScatterCreateToAll(currPast,&ctx,&(currPastSeq));
  VecScatterBegin(ctx,currPast,currPastSeq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,currPast,currPastSeq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);

  return ierr;

}
//==============================================================================

//==============================================================================
/// Assert steady state boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolverBase::assertNBC(int iR,int iZ,int iEq)
{
  if (reflectingBCs)
    assertNCurrentBC(iR,iZ,iEq);
  else if (goldinBCs)
    assertNGoldinBC(iR,iZ,iEq);
  else if (diffusionBCs)
    assertNGoldinP1BC(iR,iZ,iEq);
  else
    assertNFluxBC(iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert steady state boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolverBase::assertSBC(int iR,int iZ,int iEq)
{
  if (reflectingBCs)
    assertSCurrentBC(iR,iZ,iEq);
  else if (goldinBCs)
    assertSGoldinBC(iR,iZ,iEq);
  else if (diffusionBCs)
    assertSGoldinP1BC(iR,iZ,iEq);
  else
    assertSFluxBC(iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert steady state boundary condition on the west face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolverBase::assertWBC(int iR,int iZ,int iEq)
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
/// Assert steady state boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolverBase::assertEBC(int iR,int iZ,int iEq)
{
  if (reflectingBCs)
    assertECurrentBC(iR,iZ,iEq);
  else if (goldinBCs)
    assertEGoldinBC(iR,iZ,iEq);
  else if (diffusionBCs)
    assertEGoldinP1BC(iR,iZ,iEq);
  else
    assertEFluxBC(iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert the steady state current boundary condition on the north face at 
/// location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolverBase::assertNCurrentBC(int iR,int iZ,int iEq)
{
  northCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert the steady state current boundary condition on the south face at 
/// location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolverBase::assertSCurrentBC(int iR,int iZ,int iEq)
{
  southCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert the steady state current boundary condition on the west face at 
/// location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolverBase::assertWCurrentBC(int iR,int iZ,int iEq)
{
  westCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert the steady state current boundary condition on the south face at 
/// location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolverBase::assertECurrentBC(int iR,int iZ,int iEq)
{
  eastCurrent(1,iR,iZ,iEq);
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the north face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolverBase::assertNGoldinBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;

  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->nOutwardCurrToFluxRatioBC(iR);
  double inFluxWeightRatio = GGQD->nOutwardCurrToFluxRatioInwardWeightedBC(iR);
  double absCurrent = GGQD->nAbsCurrentBC(iR);
  double inwardCurrent = GGQD->nInwardCurrentBC(iR);
  double inwardFlux = GGQD->nInwardFluxBC(iR);

  northCurrent(1.0,iR,iZ,iEq);
  value = -ratio; 
  ierr = MatSetValue(*A,iEq,indices[iNF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = inwardCurrent-inFluxWeightRatio*inwardFlux; 
  ierr = VecSetValue(*b,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;

};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the south face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolverBase::assertSGoldinBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;

  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->sOutwardCurrToFluxRatioBC(iR);
  double inFluxWeightRatio = GGQD->sOutwardCurrToFluxRatioInwardWeightedBC(iR);
  double absCurrent = GGQD->sAbsCurrentBC(iR);
  double inwardCurrent = GGQD->sInwardCurrentBC(iR);
  double inwardFlux = GGQD->sInwardFluxBC(iR);

  southCurrent(1.0,iR,iZ,iEq);
  value = -ratio; 
  ierr = MatSetValue(*A,iEq,indices[iSF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = inwardCurrent-inFluxWeightRatio*inwardFlux; 
  ierr = VecSetValue(*b,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;

};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the east face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolverBase::assertEGoldinBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;

  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->eOutwardCurrToFluxRatioBC(iZ);
  double inFluxWeightRatio = GGQD->eOutwardCurrToFluxRatioInwardWeightedBC(iZ);
  double absCurrent = GGQD->eAbsCurrentBC(iZ);
  double inwardCurrent = GGQD->eInwardCurrentBC(iZ);
  double inwardFlux = GGQD->eInwardFluxBC(iZ);

  eastCurrent(1.0,iR,iZ,iEq);
  value = -ratio; 
  ierr = MatSetValue(*A,iEq,indices[iEF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = inwardCurrent-inFluxWeightRatio*inwardFlux; 
  ierr = VecSetValue(*b,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;

};
//==============================================================================

//==============================================================================
/// Assert Gol'din's P1 boundary condition on the north face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolverBase::assertNGoldinP1BC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  northCurrent(1.0,iR,iZ,iEq);
  value = 1.0/sqrt(3.0); 
  ierr = MatSetValue(*A,iEq,indices[iNF],value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;

};
//==============================================================================

//==============================================================================
/// Assert Gol'din's P1 boundary condition on the south face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolverBase::assertSGoldinP1BC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  southCurrent(1.0,iR,iZ,iEq);
  value = -1.0/sqrt(3.0); 
  ierr = MatSetValue(*A,iEq,indices[iSF],value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;

};
//==============================================================================

//==============================================================================
/// Assert Gol'din's P1 boundary condition on the east face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolverBase::assertEGoldinP1BC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  eastCurrent(1.0,iR,iZ,iEq);
  value = -1.0/sqrt(3.0); 
  ierr = MatSetValue(*A,iEq,indices[iSF],value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;

};
//==============================================================================

//==============================================================================
/// Calculate multigroup source coefficient for cell at (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
double GreyGroupSolverBase::calcScatterAndFissionCoeff(int iR,int iZ)
{

  double localFissionSource,localSigF,localNu,localChiP,localSigS,\
    sourceCoefficient;

  localSigS = materials->oneGroupXS->sigS(iZ,iR);
  localFissionSource = materials->oneGroupXS->qdFluxCoeff(iZ,iR);
  sourceCoefficient = localSigS + localFissionSource;

  return sourceCoefficient;
};
//==============================================================================

//==============================================================================
/// Form a portion of the linear system that belongs to GGQD 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] rEval radial location to evaluate integrating factor at
double GreyGroupSolverBase::calcIntegratingFactor(int iR,int iZ,double rEval,int iLoc) 
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
vector<int> GreyGroupSolverBase::getIndices(int iR,int iZ)
{

  vector<int> oneGroupIndices;

  // Get indices for a single energy group 
  oneGroupIndices = mesh->getQDCellIndices(iR,iZ);

  vector<int> indices {oneGroupIndices[iCF] + GGQD->indexOffset, 
                       oneGroupIndices[iWF] + GGQD->indexOffset, 
                       oneGroupIndices[iEF] + GGQD->indexOffset, 
                       oneGroupIndices[iNF] + GGQD->indexOffset, 
                       oneGroupIndices[iSF] + GGQD->indexOffset,
                       oneGroupIndices[iWC],
                       oneGroupIndices[iEC],
                       oneGroupIndices[iNC],
                       oneGroupIndices[iSC]};

  return indices;
};

//==============================================================================

//==============================================================================
/// Return geometry parameters for the cell located at (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] gParams vector containing volume and surfaces areas of the 
///   west, east, north, and south faces, in that order.
vector<double> GreyGroupSolverBase::calcGeoParams(int iR,int iZ)
{
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

  vector<double> gParams {volume,
                          wFaceSA,
                          eFaceSA,
                          nFaceSA,
                          sFaceSA};

  return gParams;         
};
//==============================================================================

//==============================================================================
/// Return volume-averaged radial coordinate for cell with boundaries rDown and
///   rUp
/// @param [in] rDown location of left cell edge
/// @param [in] rUp location of right cell edge
/// @param [out] volAvgR volume-averaged radial coordinate 
double GreyGroupSolverBase::calcVolAvgR(double rDown,double rUp)
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
double GreyGroupSolverBase::getWestErr(int iZ,int iR)
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
double GreyGroupSolverBase::getWestErz(int iZ,int iR)
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
double GreyGroupSolverBase::getEastErr(int iZ,int iR)
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
double GreyGroupSolverBase::getEastErz(int iZ,int iR)
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
double GreyGroupSolverBase::getNorthEzz(int iZ,int iR)
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
double GreyGroupSolverBase::getNorthErz(int iZ,int iR)
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
double GreyGroupSolverBase::getSouthEzz(int iZ,int iR)
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
double GreyGroupSolverBase::getSouthErz(int iZ,int iR)
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
int GreyGroupSolverBase::getFlux()
{
  vector<int> indices;
  PetscErrorCode ierr;
  PetscScalar value[5]; 
  PetscInt index[5]; 
  VecScatter     ctx;
  Vec temp_x_seq;

  // Gather values of x on all procs
  VecScatterCreateToAll(*x,&ctx,&temp_x_seq);
  VecScatterBegin(ctx,*x,temp_x_seq,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,*x,temp_x_seq,INSERT_VALUES,SCATTER_FORWARD);

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      indices = getIndices(iR,iZ);

      index[0] = indices[iCF];
      index[1] = indices[iWF];
      index[2] = indices[iEF];
      index[3] = indices[iNF];
      index[4] = indices[iSF];

      // Read fluxes into flux vector
      ierr = VecGetValues(temp_x_seq,5,index,value);
      ierr = VecSetValues(xFlux,5,index,value,INSERT_VALUES);CHKERRQ(ierr); 

      // Read fluxes into GGQD object
      GGQD->sFlux(iZ,iR) = value[0]; 
      GGQD->sFluxR(iZ,iR) = value[1];
      GGQD->sFluxR(iZ,iR+1) = value[2];
      GGQD->sFluxZ(iZ,iR) = value[3]; 
      GGQD->sFluxZ(iZ+1,iR) = value[4]; 

    }
  } 

  /* Destroy scatter context */
  VecScatterDestroy(&ctx);
  VecDestroy(&temp_x_seq);

  /* Finalize assembly of xFlux */
  VecAssemblyBegin(xFlux);
  VecAssemblyEnd(xFlux);

  return ierr;

};
//==============================================================================

//==============================================================================
/// Extract cell average values from solution vector and store
int GreyGroupSolverBase::getCurrent()
{
  vector<int> indices;
  PetscErrorCode ierr;
  PetscScalar value[4]; 
  PetscInt index[4]; 
  VecScatter     ctx;
  Vec temp_currPast_seq;

  // Gather values of x on all procs
  VecScatterCreateToAll(currPast,&ctx,&temp_currPast_seq);
  VecScatterBegin(ctx,currPast,temp_currPast_seq,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,currPast,temp_currPast_seq,INSERT_VALUES,SCATTER_FORWARD);

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      indices = getIndices(iR,iZ);

      index[0] = indices[iWC];
      index[1] = indices[iEC];
      index[2] = indices[iNC];
      index[3] = indices[iSC];

      // Read fluxes into flux vector
      ierr = VecGetValues(temp_currPast_seq,4,index,value);

      // Get cell currents 
      GGQD->currentR(iZ,iR) = value[0];
      GGQD->currentR(iZ,iR+1) = value[1]; 
      GGQD->currentZ(iZ,iR) = value[2];
      GGQD->currentZ(iZ+1,iR) = value[3];

    }
  } 

  /* Destroy scatter context */
  VecScatterDestroy(&ctx);
  VecDestroy(&temp_currPast_seq);

  return ierr;

};
//==============================================================================

//==============================================================================
/// Map values from 2D matrix to 1D solution vector
int GreyGroupSolverBase::setFlux()
{
  vector<int> indices;
  PetscErrorCode ierr;
  PetscScalar value[5]; 
  PetscInt index[5]; 
  VecScatter     ctx;

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      indices = getIndices(iR,iZ);

      value[0] = GGQD->sFlux(iZ,iR);
      value[1] = GGQD->sFluxR(iZ,iR);
      value[2] = GGQD->sFluxR(iZ,iR+1);
      value[3] = GGQD->sFluxZ(iZ,iR);
      value[4] = GGQD->sFluxZ(iZ+1,iR);

      index[0] = indices[iCF];
      index[1] = indices[iWF]; 
      index[2] = indices[iEF]; 
      index[3] = indices[iNF]; 
      index[4] = indices[iSF]; 

      ierr = VecSetValues(&xPast,5,index,value,INSERT_VALUES);CHKERRQ(ierr); 
    }
  }  

  // Finalize assembly of xPast
  VecAssemblyBegin(*xPast);
  VecAssemblyEnd(*xPast);

  // Form xPast_seq
  VecScatterCreateToAll(*xPast,&ctx,xPastSeq);
  VecScatterBegin(ctx,*xPast,*xPastSeq,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,*xPast,*xPastSeq,INSERT_VALUES,SCATTER_FORWARD);

  /* Destroy scatter context */
  VecScatterDestroy(&ctx);


  return ierr;

};
//==============================================================================

//==============================================================================

/// Extract flux values from GGQD object, store to solution vector, and return 
/// @param [out] solVector flux values mapped to the 1D vector
Eigen::VectorXd GreyGroupSolverBase::getFluxSolutionVector()
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
Eigen::VectorXd GreyGroupSolverBase::getCurrentSolutionVector()
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
/// Assert the flux boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolverBase::assertNFluxBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  value = 1.0;
  ierr = MatSetValue(*A,iEq,indices[iNF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = GGQD->nFluxBC(iR);
  ierr = VecSetValue(*b,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  //Atemp.insert(iEq,indices[iNF]) = 1.0;
  //(*b)(iEq) = GGQD->nFluxBC(iR);

  return ierr;

};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolverBase::assertSFluxBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  value = 1.0;
  ierr = MatSetValue(*A,iEq,indices[iSF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = GGQD->sFluxBC(iR);
  ierr = VecSetValue(*b,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  //Atemp.insert(iEq,indices[iSF]) = 1.0;
  //(*b)(iEq) = GGQD->sFluxBC(iR);

  return ierr;

};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the west face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolverBase::assertWFluxBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  value = 1.0;
  ierr = MatSetValue(*A,iEq,indices[iWF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = GGQD->wFluxBC(iZ);
  ierr = VecSetValue(*b,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  //Atemp.insert(iEq,indices[iWF]) = 1.0;
  //(*b)(iEq) = GGQD->wFluxBC(iZ);

  return ierr;

};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolverBase::assertEFluxBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  value = 1.0;
  ierr = MatSetValue(*A,iEq,indices[iEF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = GGQD->eFluxBC(iZ);
  ierr = VecSetValue(*b,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  //Atemp.insert(iEq,indices[iEF]) = 1.0;
  //(*b)(iEq) = GGQD->eFluxBC(iZ);

  return ierr;

};
//==============================================================================

//==============================================================================
/// Assign pointers to linear system components 
///
void GreyGroupSolverBase::assignMPQDPointer(MultiPhysicsCoupledQD * myMPQD)
{

  MPQD = myMPQD;

};
//==============================================================================

//==============================================================================
/// Assign pointers to linear system components 
///
void GreyGroupSolverBase::assignLinearSystemPointers(Mat * myA,\
    Vec * myx,\
    Vec * myxPast,\
    Vec * myxPastSeq,\
    Vec * myb)
{

  A = myA;
  x = myx;
  xPast = myxPast;
  xPastSeq = myxPastSeq;
  b = myb;

};
//==============================================================================


//==============================================================================
/// Check for optional inputs of relevance to this object
void GreyGroupSolverBase::checkOptionalParams()
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
    else if (boundaryType == "diffusion" or boundaryType == "DIFFUSION" \
        or boundaryType == "Diffusion")
      diffusionBCs = true;


  }
}
//==============================================================================
