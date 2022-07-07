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

  // initialize size of linear system
  C.resize(nCurrentUnknowns,nUnknowns);
  C.reserve(4*nCurrentUnknowns);
  xFlux.setZero(nUnknowns);
  currPast.setZero(nCurrentUnknowns);
  d.setZero(nCurrentUnknowns);

  /* Initialize PETSc variables */
  // Multiphysics system variables 
  initPETScRectMat(&C_p,nCurrentUnknowns,nUnknowns,20);
  initPETScVec(&currPast_p,nCurrentUnknowns);
  initPETScVec(&d_p,nCurrentUnknowns);
  initPETScVec(&xFlux_p,nUnknowns);

  // Broadcast currPast
  VecScatterCreateToAll(currPast_p,&ctx,&(currPast_p_seq));
  VecScatterBegin(ctx,currPast_p,currPast_p_seq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,currPast_p,currPast_p_seq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);

  // Set number of unknowns if GreyGroupQD object
  GGQD->nUnknowns = nUnknowns;
  GGQD->nCurrentUnknowns = nCurrentUnknowns;

  checkOptionalParams();
};

//==============================================================================

//==============================================================================
/// Calculate multigroup source coefficient for cell at (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
double GreyGroupSolver::calcScatterAndFissionCoeff(int iR,int iZ)
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
int GreyGroupSolver::getFlux()
{
  vector<int> indices;
  PetscErrorCode ierr;
  PetscScalar value[5]; 
  PetscInt index[5]; 
  VecScatter     ctx;
  Vec temp_x_p_seq;

  // Gather values of x_p on all procs
  VecScatterCreateToAll(MPQD->x_p,&ctx,&temp_x_p_seq);
  VecScatterBegin(ctx,MPQD->x_p,temp_x_p_seq,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,MPQD->x_p,temp_x_p_seq,INSERT_VALUES,SCATTER_FORWARD);

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
      ierr = VecGetValues(temp_x_p_seq,5,index,value);
      ierr = VecSetValues(xFlux_p,5,index,value,INSERT_VALUES);CHKERRQ(ierr); 

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
  VecDestroy(&temp_x_p_seq);

  /* Finalize assembly of xFlux_p */
  VecAssemblyBegin(xFlux_p);
  VecAssemblyEnd(xFlux_p);

  return ierr;
};
//==============================================================================

//==============================================================================
/// Extract cell average values from solution vector and store
int GreyGroupSolver::getCurrent()
{
  vector<int> indices;
  PetscErrorCode ierr;
  PetscScalar value[4]; 
  PetscInt index[4]; 
  VecScatter     ctx;
  Vec temp_currPast_p_seq;

  // Gather values of x_p on all procs
  VecScatterCreateToAll(currPast_p,&ctx,&temp_currPast_p_seq);
  VecScatterBegin(ctx,currPast_p,temp_currPast_p_seq,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,currPast_p,temp_currPast_p_seq,INSERT_VALUES,SCATTER_FORWARD);

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
      ierr = VecGetValues(temp_currPast_p_seq,4,index,value);

      // Get cell currents 
      GGQD->currentR(iZ,iR) = value[0];
      GGQD->currentR(iZ,iR+1) = value[1]; 
      GGQD->currentZ(iZ,iR) = value[2];
      GGQD->currentZ(iZ+1,iR) = value[3];

    }
  } 

  /* Destroy scatter context */
  VecScatterDestroy(&ctx);
  VecDestroy(&temp_currPast_p_seq);

  return ierr;
};
//==============================================================================

//==============================================================================
/// Map values from 2D matrix to 1D solution vector
int GreyGroupSolver::setFlux()
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

      ierr = VecSetValues(MPQD->xPast_p,5,index,value,INSERT_VALUES);CHKERRQ(ierr); 
    }
  }  

  // Finalize assembly of xPast_p
  VecAssemblyBegin(MPQD->xPast_p);
  VecAssemblyEnd(MPQD->xPast_p);

  // Form xPast_p_seq
  VecScatterCreateToAll(MPQD->xPast_p,&ctx,&(MPQD->xPast_p_seq));
  VecScatterBegin(ctx,MPQD->xPast_p,MPQD->xPast_p_seq,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,MPQD->xPast_p,MPQD->xPast_p_seq,INSERT_VALUES,SCATTER_FORWARD);

  /* Destroy scatter context */
  VecScatterDestroy(&ctx);


  return ierr;
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

/* PETSc functions */

/* STEADY STATE */

//==============================================================================
/// Form a portion of the linear system  
void GreyGroupSolver::formSteadyStateLinearSystem()	      
{

  int iEq = GGQD->indexOffset;

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
};
//==============================================================================

//==============================================================================
/// Form a portion of the current back calc linear system that belongs to GGQD 
///
int GreyGroupSolver::formSteadyStateBackCalcSystem()	      
{
  int iEq = GGQD->indexOffset;
  PetscErrorCode ierr;

  // Reset linear system
  MatZeroEntries(C_p);
  VecZeroEntries(d_p);

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

  /* Finalize assembly for C_p and d_p */
  ierr = MatAssemblyBegin(C_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  
  ierr = VecAssemblyBegin(d_p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(d_p);CHKERRQ(ierr);

  return ierr;

};

//==============================================================================


//==============================================================================
/// Assert the steady state zeroth moment equation for cell (iR,iZ)
///
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertSteadyStateZerothMoment(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double vPast = materials->oneGroupXS->neutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double scatterCoeff, fissionCoeff, keff, cellFlux;
  PetscErrorCode ierr;
  PetscScalar value;
  PetscInt index;

  indices = getIndices(iR,iZ);

  // Scattering source term (implicit)
  scatterCoeff = materials->oneGroupXS->sigS(iZ,iR);
  value = -geoParams[iCF] * scatterCoeff; 
  index = indices[iCF];
  ierr = MatSetValue(MPQD->A_p,iEq,index,value,ADD_VALUES);CHKERRQ(ierr); 

  // Fission source term (explicit for power iteration)
  fissionCoeff = materials->oneGroupXS->qdFluxCoeff(iZ,iR);
  keff = materials->oneGroupXS->keff;
  cellFlux = GGQD->sFlux(iZ,iR);
  value = geoParams[iCF]*(fissionCoeff*cellFlux/keff + GGQD->q(iZ,iR));
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  // DNP source term
  GGQD->mpqd->dnpSource(iZ,iR,iEq,-geoParams[iCF]);

  // populate entries representing streaming and reaction terms
  value = geoParams[iCF] * sigT;
  index = indices[iCF];
  ierr = MatSetValue(MPQD->A_p,iEq,index,value,ADD_VALUES);CHKERRQ(ierr); 

  steadyStateWestCurrent(-geoParams[iWF],iR,iZ,iEq);

  steadyStateEastCurrent(geoParams[iEF],iR,iZ,iEq);

  steadyStateNorthCurrent(-geoParams[iNF],iR,iZ,iEq);

  steadyStateSouthCurrent(geoParams[iSF],iR,iZ,iEq);

  return ierr;
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
/// Enforce coefficients for steady state current on south face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::steadyStateSouthCurrent(double coeff,int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];

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

  index[0] = indices[iSF]; value[0] = -coeff*(EzzS/deltaZ+zetaL);

  index[1] = indices[iCF]; value[1] = coeff*EzzC/deltaZ;

  index[2] = indices[iWF]; value[2] = coeff*(rDown*ErzW/(rAvg*deltaR));

  index[3] = indices[iEF]; value[3] = -coeff*rUp*ErzE/(rAvg*deltaR);

  ierr = MatSetValues(MPQD->A_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  return ierr;
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for steady state current on north face 
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::steadyStateNorthCurrent(double coeff,int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];

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

  index[0] = indices[iNF]; value[0] = coeff*(EzzN/deltaZ - zetaL);

  index[1] = indices[iCF]; value[1] = -coeff*EzzC/deltaZ;

  index[2] = indices[iWF]; value[2] = coeff*(rDown*ErzW/(rAvg*deltaR));

  index[3] = indices[iEF]; value[3] = -coeff*rUp*ErzE/(rAvg*deltaR);

  ierr = MatSetValues(MPQD->A_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  return ierr;
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for steady state current on west face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::steadyStateWestCurrent(double coeff,int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];

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

  index[0] = indices[iSF]; value[0] = -coeff*ErzS/deltaZ;

  index[1] = indices[iNF]; value[1] = coeff*ErzN/deltaZ;

  index[2] = indices[iCF]; value[2] = -coeff*hCent*ErrC/(hDown*deltaR);

  index[3] = indices[iWF]; value[3] = coeff*(hDown*ErrW/(hDown*deltaR) - zetaL);

  ierr = MatSetValues(MPQD->A_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  return ierr;
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for steady state current on east face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::steadyStateEastCurrent(double coeff,int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];

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

  index[0] = indices[iSF]; value[0] = -coeff*ErzS/deltaZ;

  index[1] = indices[iNF]; value[1] = coeff*ErzN/deltaZ;

  index[2] = indices[iCF]; value[2] = coeff*hCent*ErrC/(hUp*deltaR);

  index[3] = indices[iEF]; value[3] = -coeff*(hUp*ErrE/(hUp*deltaR) + zetaL);

  ierr = MatSetValues(MPQD->A_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  return ierr;
};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate steady state current on south face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::calcSteadyStateSouthCurrent(int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];


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

  index[0] = indices[iSF]; value[0] = -coeff*(EzzS/deltaZ + zetaL);

  index[1] = indices[iCF]; value[1] = coeff*EzzC/deltaZ;

  index[2] = indices[iWF]; value[2] = coeff*(rDown*ErzW/(rAvg*deltaR));

  index[3] = indices[iEF]; value[3] = -coeff*rUp*ErzE/(rAvg*deltaR);

  ierr = MatSetValues(C_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  return ierr;
};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate steady state current on north face 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::calcSteadyStateNorthCurrent(int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];


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

  index[0] = indices[iNF]; value[0] = coeff*(EzzN/deltaZ - zetaL);

  index[1] = indices[iCF]; value[1] = -coeff*EzzC/deltaZ;

  index[2] = indices[iWF]; value[2] = coeff*(rDown*ErzW/(rAvg*deltaR));

  index[3] = indices[iEF]; value[3] = -coeff*rUp*ErzE/(rAvg*deltaR);

  ierr = MatSetValues(C_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  return ierr;
};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate steady state current on west face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::calcSteadyStateWestCurrent(int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];

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

  index[0] = indices[iSF]; value[0] = -coeff*ErzS/deltaZ;

  index[1] = indices[iNF]; value[1] = coeff*ErzN/deltaZ;

  index[2] = indices[iCF]; value[2] = -coeff*hCent*ErrC/(hDown*deltaR);

  index[3] = indices[iWF]; value[3] = coeff*(hDown*ErrW/(hDown*deltaR) - zetaL);

  ierr = MatSetValues(C_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  return ierr;
};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate steady state current on east face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::calcSteadyStateEastCurrent(int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];

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

  index[0] = indices[iSF]; value[0] = -coeff*ErzS/deltaZ;

  index[1] = indices[iNF]; value[1] = coeff*ErzN/deltaZ;

  index[2] = indices[iCF]; value[2] = coeff*hCent*ErrC/(hUp*deltaR);

  index[3] = indices[iEF]; value[3] = -coeff*(hUp*ErrE/(hUp*deltaR) + zetaL);

  ierr = MatSetValues(C_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  return ierr;
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
  else if (diffusionBCs)
    assertSteadyStateNGoldinP1BC(iR,iZ,iEq);
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
  else if (diffusionBCs)
    assertSteadyStateSGoldinP1BC(iR,iZ,iEq);
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
  else if (diffusionBCs)
    assertSteadyStateEGoldinP1BC(iR,iZ,iEq);
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
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertSteadyStateNGoldinBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;

  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->nOutwardCurrToFluxRatioBC(iR);
  double inFluxWeightRatio = GGQD->nOutwardCurrToFluxRatioInwardWeightedBC(iR);
  double absCurrent = GGQD->nAbsCurrentBC(iR);
  double inwardCurrent = GGQD->nInwardCurrentBC(iR);
  double inwardFlux = GGQD->nInwardFluxBC(iR);

  steadyStateNorthCurrent(1.0,iR,iZ,iEq);
  value = -ratio; 
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iNF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = inwardCurrent-inFluxWeightRatio*inwardFlux; 
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the south face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertSteadyStateSGoldinBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;

  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->sOutwardCurrToFluxRatioBC(iR);
  double inFluxWeightRatio = GGQD->sOutwardCurrToFluxRatioInwardWeightedBC(iR);
  double absCurrent = GGQD->sAbsCurrentBC(iR);
  double inwardCurrent = GGQD->sInwardCurrentBC(iR);
  double inwardFlux = GGQD->sInwardFluxBC(iR);

  steadyStateSouthCurrent(1.0,iR,iZ,iEq);
  value = -ratio; 
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iSF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = inwardCurrent-inFluxWeightRatio*inwardFlux; 
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the east face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertSteadyStateEGoldinBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;

  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->eOutwardCurrToFluxRatioBC(iZ);
  double inFluxWeightRatio = GGQD->eOutwardCurrToFluxRatioInwardWeightedBC(iZ);
  double absCurrent = GGQD->eAbsCurrentBC(iZ);
  double inwardCurrent = GGQD->eInwardCurrentBC(iZ);
  double inwardFlux = GGQD->eInwardFluxBC(iZ);

  steadyStateEastCurrent(1.0,iR,iZ,iEq);
  value = -ratio; 
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iEF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = inwardCurrent-inFluxWeightRatio*inwardFlux; 
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's P1 boundary condition on the north face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertSteadyStateNGoldinP1BC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  steadyStateNorthCurrent(1.0,iR,iZ,iEq);
  value = 1.0 / sqrt(3.0);
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iNF],value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's P1 boundary condition on the south face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertSteadyStateSGoldinP1BC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  steadyStateSouthCurrent(1.0,iR,iZ,iEq);
  value = -1.0 / sqrt(3.0);
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iSF],value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's P1 boundary condition on the east face at location (iR,iZ)
/// for steady state solves
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertSteadyStateEGoldinP1BC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  steadyStateEastCurrent(1.0,iR,iZ,iEq);
  value = -1.0 / sqrt(3.0); 
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iEF],value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

/* TRANSIENT AND STEADY STATE */

//==============================================================================
/// Assert the flux boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertNFluxBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  value = 1.0;
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iNF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = GGQD->nFluxBC(iR);
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertSFluxBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  value = 1.0;
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iSF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = GGQD->sFluxBC(iR);
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the west face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertWFluxBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  value = 1.0;
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iWF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = GGQD->wFluxBC(iZ);
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertEFluxBC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  value = 1.0;
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iEF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = GGQD->eFluxBC(iZ);
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Compute currents from flux values in x
///
int GreyGroupSolver::backCalculateCurrent()
{
  // PETSc object to broadcast variables 
  VecScatter     ctx;
  PetscErrorCode ierr;

  ierr = MatMultAdd(C_p,xFlux_p,d_p,currPast_p);

  // Broadcast currPast
  VecDestroy(&(currPast_p_seq));
  VecScatterCreateToAll(currPast_p,&ctx,&(currPast_p_seq));
  VecScatterBegin(ctx,currPast_p,currPast_p_seq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,currPast_p,currPast_p_seq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);

  return ierr;
}
//==============================================================================

/* TRANSIENT */

//==============================================================================
/// Form a portion of the linear system  

void GreyGroupSolver::formLinearSystem()	      
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

//==============================================================================
/// Assert the transient zeroth moment equation for cell (iR,iZ)
///
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertZerothMoment(int iR,int iZ,int iEq)
{
  vector<int> indices;
  vector<double> geoParams = mesh->getGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->oneGroupXS->neutV(iZ,iR);
  double vPast = materials->oneGroupXS->neutVPast(iZ,iR);
  double sigT = materials->oneGroupXS->sigT(iZ,iR);
  double groupSourceCoeff;
  PetscErrorCode ierr;
  PetscScalar value,past_flux;
  PetscInt index;

  indices = getIndices(iR,iZ);

  // Scatter and fission source term
  groupSourceCoeff = calcScatterAndFissionCoeff(iR,iZ);
  value = -geoParams[iCF] * groupSourceCoeff; 
  index = indices[iCF];
  ierr = MatSetValue(MPQD->A_p,iEq,index,value,ADD_VALUES);CHKERRQ(ierr); 

  // DNP source term
  GGQD->mpqd->dnpSource(iZ,iR,iEq,-geoParams[iCF]);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  value = geoParams[iCF] * ((1/(v*deltaT)) + sigT);
  index = indices[iCF];
  ierr = MatSetValue(MPQD->A_p,iEq,index,value,ADD_VALUES);CHKERRQ(ierr); 

  westCurrent(-geoParams[iWF],iR,iZ,iEq);

  eastCurrent(geoParams[iEF],iR,iZ,iEq);

  northCurrent(-geoParams[iNF],iR,iZ,iEq);

  southCurrent(geoParams[iSF],iR,iZ,iEq);

  // formulate RHS entry
  VecGetValues(MPQD->xPast_p_seq,1,&index,&past_flux);CHKERRQ(ierr);
  value = geoParams[iCF]*((past_flux/(vPast*deltaT)) + GGQD->q(iZ,iR));
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
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
/// Enforce coefficients for transient current on south face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::southCurrent(double coeff,int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];
  double curr_value;
  int curr_index;

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

  index[0] = indices[iSF]; value[0] = -coeff*(EzzS/deltaZ+zetaL);

  index[1] = indices[iCF]; value[1] = coeff*EzzC/deltaZ;

  index[2] = indices[iWF]; value[2] = coeff*(rDown*ErzW/(rAvg*deltaR));

  index[3] = indices[iEF]; value[3] = -coeff*rUp*ErzE/(rAvg*deltaR);

  ierr = MatSetValues(MPQD->A_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {

      // Get index of current
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      curr_index = indices[iSC];

      // Get previous multigroup current and neutron velocity
      VecGetValues(GGQD->mgqd->QDSolve->currPast_p_seq,1,&curr_index,&curr_value);CHKERRQ(ierr);
      mgqdNeutV = materials->zNeutVel(iZ+1,iR,iGroup); 

      // Set coefficient (curr_value) in RHS vector
      curr_value = -(coeff/deltaT)*(curr_value/mgqdNeutV);;
      ierr = VecSetValue(MPQD->b_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
    }
  }
  else
  {
    curr_index = indices[iSC];
    VecGetValues(currPast_p_seq,1,&curr_index,&curr_value);CHKERRQ(ierr);
    curr_value = -coeff*(curr_value/(vPast*deltaT));
    ierr = VecSetValue(MPQD->b_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
  }

  return ierr;
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for transient current on north face 
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::northCurrent(double coeff,int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];
  double curr_value;
  int curr_index;

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

  index[0] = indices[iNF]; value[0] = coeff*(EzzN/deltaZ - zetaL);

  index[1] = indices[iCF]; value[1] = -coeff*EzzC/deltaZ;

  index[2] = indices[iWF]; value[2] = coeff*(rDown*ErzW/(rAvg*deltaR));

  index[3] = indices[iEF]; value[3] = -coeff*rUp*ErzE/(rAvg*deltaR);

  ierr = MatSetValues(MPQD->A_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      // Get index of current
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      curr_index = indices[iNC];

      // Get previous multigroup current and neutron velocity
      VecGetValues(GGQD->mgqd->QDSolve->currPast_p_seq,1,&curr_index,&curr_value);CHKERRQ(ierr);
      mgqdNeutV = materials->zNeutVel(iZ,iR,iGroup); 

      // Set coefficient (curr_value) in RHS vector
      curr_value = -(coeff/deltaT)*(curr_value/mgqdNeutV);;
      ierr = VecSetValue(MPQD->b_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
    }
  }
  else
  {
    curr_index = indices[iNC];
    VecGetValues(currPast_p_seq,1,&curr_index,&curr_value);CHKERRQ(ierr);
    curr_value = -coeff*(curr_value/(vPast*deltaT));
    ierr = VecSetValue(MPQD->b_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
  }

  return ierr;

};
//==============================================================================

//==============================================================================
/// Enforce coefficients for transient current on west face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::westCurrent(double coeff,int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];
  double curr_value;
  int curr_index;

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

  index[0] = indices[iSF]; value[0] = -coeff*ErzS/deltaZ;

  index[1] = indices[iNF]; value[1] = coeff*ErzN/deltaZ;

  index[2] = indices[iCF]; value[2] = -coeff*hCent*ErrC/(hDown*deltaR);

  index[3] = indices[iWF]; value[3] = coeff*(hDown*ErrW/(hDown*deltaR) - zetaL);

  ierr = MatSetValues(MPQD->A_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      // Get index of current
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      curr_index = indices[iWC];

      // Get previous multigroup current and neutron velocity
      VecGetValues(GGQD->mgqd->QDSolve->currPast_p_seq,1,&curr_index,&curr_value);CHKERRQ(ierr);
      mgqdNeutV = materials->rNeutVel(iZ,iR,iGroup); 

      // Set coefficient (curr_value) in RHS vector
      curr_value = -(coeff/deltaT)*(curr_value/mgqdNeutV);
      ierr = VecSetValue(MPQD->b_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
    }
  }
  else
  {
    curr_index = indices[iWC];
    VecGetValues(currPast_p_seq,1,&curr_index,&curr_value);CHKERRQ(ierr);
    curr_value = -coeff*(curr_value/(vPast*deltaT));
    ierr = VecSetValue(MPQD->b_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
  }

  return ierr;
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for transient current on east face
///
/// @param [in] coeff coefficient multiplying current
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::eastCurrent(double coeff,int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar value[4];
  PetscInt index[4];
  double curr_value;
  int curr_index;

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

  index[0] = indices[iSF]; value[0] = -coeff*ErzS/deltaZ;

  index[1] = indices[iNF]; value[1] = coeff*ErzN/deltaZ;

  index[2] = indices[iCF]; value[2] = coeff*hCent*ErrC/(hUp*deltaR);

  index[3] = indices[iEF]; value[3] = -coeff*(hUp*ErrE/(hUp*deltaR) + zetaL);

  ierr = MatSetValues(MPQD->A_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      // Get index of current
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      curr_index = indices[iEC];

      // Get previous multigroup current and neutron velocity
      VecGetValues(GGQD->mgqd->QDSolve->currPast_p_seq,1,&curr_index,&curr_value);CHKERRQ(ierr);
      mgqdNeutV = materials->rNeutVel(iZ,iR+1,iGroup); 

      // Set coefficient (curr_value) in RHS vector
      curr_value = -(coeff/deltaT)*(curr_value/mgqdNeutV);;
      ierr = VecSetValue(MPQD->b_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
    }
  }
  else
  {
    curr_index = indices[iEC];
    VecGetValues(currPast_p_seq,1,&curr_index,&curr_value);CHKERRQ(ierr);
    curr_value = -coeff*(curr_value/(vPast*deltaT));
    ierr = VecSetValue(MPQD->b_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 

  }

  return ierr;
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
void GreyGroupSolver::assertSBC(int iR,int iZ,int iEq)
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
/// Assert steady state boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
void GreyGroupSolver::assertEBC(int iR,int iZ,int iEq)
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
int GreyGroupSolver::assertNGoldinBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->nOutwardCurrToFluxRatioBC(iR);
  double inFluxWeightRatio = GGQD->nOutwardCurrToFluxRatioInwardWeightedBC(iR);
  double absCurrent = GGQD->nAbsCurrentBC(iR);
  double inwardCurrent = GGQD->nInwardCurrentBC(iR);
  double inwardFlux = GGQD->nInwardFluxBC(iR);
  PetscErrorCode ierr;
  double value;

  northCurrent(1.0,iR,iZ,iEq);
  value = -ratio; 
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iNF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = inwardCurrent-inFluxWeightRatio*inwardFlux; 
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertSGoldinBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->sOutwardCurrToFluxRatioBC(iR);
  double inFluxWeightRatio = GGQD->sOutwardCurrToFluxRatioInwardWeightedBC(iR);
  double absCurrent = GGQD->sAbsCurrentBC(iR);
  double inwardCurrent = GGQD->sInwardCurrentBC(iR);
  double inwardFlux = GGQD->sInwardFluxBC(iR);
  PetscErrorCode ierr;
  double value;

  southCurrent(1.0,iR,iZ,iEq);
  value = -ratio; 
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iSF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = inwardCurrent-inFluxWeightRatio*inwardFlux; 
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertEGoldinBC(int iR,int iZ,int iEq)
{
  vector<int> indices = getIndices(iR,iZ);
  double ratio = GGQD->eOutwardCurrToFluxRatioBC(iZ);
  double inFluxWeightRatio = GGQD->eOutwardCurrToFluxRatioInwardWeightedBC(iZ);
  double absCurrent = GGQD->eAbsCurrentBC(iZ);
  double inwardCurrent = GGQD->eInwardCurrentBC(iZ);
  double inwardFlux = GGQD->eInwardFluxBC(iZ);
  PetscErrorCode ierr;
  double value;

  eastCurrent(1.0,iR,iZ,iEq);
  value = -ratio; 
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iEF],value,ADD_VALUES);CHKERRQ(ierr); 
  value = inwardCurrent-inFluxWeightRatio*inwardFlux; 
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's P1 boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertNGoldinP1BC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  northCurrent(1.0,iR,iZ,iEq);
  value = 1.0/sqrt(3.0); 
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iNF],value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's P1 boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertSGoldinP1BC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  southCurrent(1.0,iR,iZ,iEq);
  value = -1.0/sqrt(3.0); 
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iSF],value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's P1 boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::assertEGoldinP1BC(int iR,int iZ,int iEq)
{
  PetscErrorCode ierr;
  double value;
  vector<int> indices = getIndices(iR,iZ);

  eastCurrent(1.0,iR,iZ,iEq);
  value = -1.0/sqrt(3.0); 
  ierr = MatSetValue(MPQD->A_p,iEq,indices[iEF],value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Form a portion of the current back calc linear system that belongs to GGQD 
///
int GreyGroupSolver::formBackCalcSystem()	      
{
  int iEq = GGQD->indexOffset;
  PetscErrorCode ierr;

  // Reset linear system
  MatZeroEntries(C_p);
  VecZeroEntries(d_p);

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

  /* Finalize assembly for C_p and d_p */
  ierr = MatAssemblyBegin(C_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(C_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  
  ierr = VecAssemblyBegin(d_p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(d_p);CHKERRQ(ierr);

  return ierr;
};

//==============================================================================


//==============================================================================
/// Enforce coefficients to calculate transient current on south face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::calcSouthCurrent(int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar curr_value,value[4];
  PetscInt curr_index,index[4];


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

  index[0] = indices[iSF]; value[0] = -coeff*(EzzS/deltaZ + zetaL);

  index[1] = indices[iCF]; value[1] = coeff*EzzC/deltaZ;

  index[2] = indices[iWF]; value[2] = coeff*(rDown*ErzW/(rAvg*deltaR));

  index[3] = indices[iEF]; value[3] = -coeff*rUp*ErzE/(rAvg*deltaR);

  ierr = MatSetValues(C_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      // Get index of current
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      curr_index = indices[iSC];

      // Get previous multigroup current and neutron velocity
      VecGetValues(GGQD->mgqd->QDSolve->currPast_p_seq,1,&curr_index,\
          &curr_value);CHKERRQ(ierr);
      mgqdNeutV = materials->zNeutVel(iZ+1,iR,iGroup); 

      // Set coefficient (curr_value) in RHS vector
      curr_value = (coeff/deltaT)*(curr_value/mgqdNeutV);;
      ierr = VecSetValue(d_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
    }
  }
  else
  {
    curr_value = coeff*(currPast(indices[iSC])/(vPast*deltaT));
    ierr = VecSetValue(d_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
  }

  return ierr;

};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate transient current on north face 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::calcNorthCurrent(int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar curr_value,value[4];
  PetscInt curr_index,index[4];

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

  index[0] = indices[iNF]; value[0] = coeff*(EzzN/deltaZ - zetaL);

  index[1] = indices[iCF]; value[1] = -coeff*EzzC/deltaZ;

  index[2] = indices[iWF]; value[2] = coeff*(rDown*ErzW/(rAvg*deltaR));

  index[3] = indices[iEF]; value[3] = -coeff*rUp*ErzE/(rAvg*deltaR);

  ierr = MatSetValues(C_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      // Get index of current
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      curr_index = indices[iNC];

      // Get previous multigroup current and neutron velocity
      VecGetValues(GGQD->mgqd->QDSolve->currPast_p_seq,\
          1,&curr_index,&curr_value);CHKERRQ(ierr);
      mgqdNeutV = materials->zNeutVel(iZ,iR,iGroup); 

      // Set coefficient (curr_value) in RHS vector
      curr_value = (coeff/deltaT)*(curr_value/mgqdNeutV);;
      ierr = VecSetValue(d_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
    }
  }
  else
  {
    curr_value = coeff*(currPast(indices[iNC])/(vPast*deltaT));
    ierr = VecSetValue(d_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
  }

  return ierr;

};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate transient current on west face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::calcWestCurrent(int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar curr_value,value[4];
  PetscInt curr_index,index[4];

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

  index[0] = indices[iSF]; value[0] = -coeff*ErzS/deltaZ;

  index[1] = indices[iNF]; value[1] = coeff*ErzN/deltaZ;

  index[2] = indices[iCF]; value[2] = -coeff*hCent*ErrC/(hDown*deltaR);

  index[3] = indices[iWF]; value[3] = coeff*(hDown*ErrW/(hDown*deltaR) - zetaL);

  ierr = MatSetValues(C_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      // Get index of current
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      curr_index = indices[iWC];

      // Get previous multigroup current and neutron velocity
      VecGetValues(GGQD->mgqd->QDSolve->currPast_p_seq,\
          1,&curr_index,&curr_value);CHKERRQ(ierr);
      mgqdNeutV = materials->rNeutVel(iZ,iR,iGroup); 

      // Set coefficient (curr_value) in RHS vector
      curr_value = (coeff/deltaT)*(curr_value/mgqdNeutV);;
      ierr = VecSetValue(d_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
    }
  }
  else
  {
    curr_value = coeff*(currPast(indices[iWC])/(vPast*deltaT));
    ierr = VecSetValue(d_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
  }

  return ierr;

};
//==============================================================================

//==============================================================================
/// Enforce coefficients to calculate transient current on east face
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
int GreyGroupSolver::calcEastCurrent(int iR,int iZ,int iEq)
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
  PetscErrorCode ierr;
  PetscScalar curr_value,value[4];
  PetscInt curr_index,index[4];

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

  index[0] = indices[iSF]; value[0] = -coeff*ErzS/deltaZ;

  index[1] = indices[iNF]; value[1] = coeff*ErzN/deltaZ;

  index[2] = indices[iCF]; value[2] = coeff*hCent*ErrC/(hUp*deltaR);

  index[3] = indices[iEF]; value[3] = -coeff*(hUp*ErrE/(hUp*deltaR) + zetaL);

  ierr = MatSetValues(C_p,1,&iEq,4,index,value,ADD_VALUES);CHKERRQ(ierr);

  // formulate RHS entry
  if (GGQD->useMGQDSources)
  {
    for (int iGroup = 0; iGroup < materials->nGroups; iGroup++)
    {
      // Get index of current
      indices = GGQD->mgqd->QDSolve->getIndices(iR,iZ,iGroup);
      curr_index = indices[iEC];

      // Get previous multigroup current and neutron velocity
      VecGetValues(GGQD->mgqd->QDSolve->currPast_p_seq,\
          1,&curr_index,&curr_value);CHKERRQ(ierr);
      mgqdNeutV = materials->rNeutVel(iZ,iR+1,iGroup); 

      // Set coefficient (curr_value) in RHS vector
      curr_value = (coeff/deltaT)*(curr_value/mgqdNeutV);;
      ierr = VecSetValue(d_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
    }
  }
  else
  {
    curr_value = coeff*(currPast(indices[iEC])/(vPast*deltaT));
    ierr = VecSetValue(d_p,iEq,curr_value,ADD_VALUES);CHKERRQ(ierr); 
  }

  return ierr;

};
//==============================================================================


//==============================================================================
/// Assign pointers to linear system components 
///
void GreyGroupSolver::assignMPQDPointer(MultiPhysicsCoupledQD * myMPQD)
{

  MPQD = myMPQD;

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
    else if (boundaryType == "diffusion" or boundaryType == "DIFFUSION" \
        or boundaryType == "Diffusion")
      diffusionBCs = true;


  }
}
//==============================================================================
