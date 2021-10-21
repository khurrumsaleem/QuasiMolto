// File: GreyGroupSolver.cpp
// Purpose: Solve RZ quasidiffusion equations  
// Date: February 05, 2020

#include "GreyGroupSolverTransient.h"
#include "GreyGroupQD.h"

using namespace std; 

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
  GGQD->mpqd->dnpSource(iZ,iR,iEq,-geoParams[iCF], &Atemp);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ);

  value = geoParams[iCF] * ((1/(v*deltaT)) + sigT);
  index = indices[iCF];
  ierr = MatSetValue(MPQD->A_p,iEq,index,value,ADD_VALUES);CHKERRQ(ierr); 

  westCurrent_p(-geoParams[iWF],iR,iZ,iEq);

  eastCurrent_p(geoParams[iEF],iR,iZ,iEq);

  northCurrent_p(-geoParams[iNF],iR,iZ,iEq);

  southCurrent_p(geoParams[iSF],iR,iZ,iEq);

  // formulate RHS entry
  VecGetValues(MPQD->xPast_p_seq,1,&index,&past_flux);CHKERRQ(ierr);
  value = geoParams[iCF]*((past_flux/(vPast*deltaT)) + GGQD->q(iZ,iR));
  ierr = VecSetValue(MPQD->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;

};
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
