// File: QuasidiffusionSolver.cpp
// Purpose: Solve RZ quasidiffusion equations  
// Date: February 05, 2020

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Materials.h"
#include "Material.h"
#include "MMS.h"
#include "QuasidiffusionSolver.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
/// QuasidiffusionSolver object constructor
///
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myInput YAML input object for the simulation
QDSolver::QDSolver(Mesh * myMesh,\
  Materials * myMaterials,\
  YAML::Node * myInput)	      
{
  // Point to variables for mesh and input file
  mesh = myMesh;
  input = myInput;
  materials = myMaterials;

  // temporary variables for initialization
  int nUnknowns;

  // calculate number of unknowns  
  energyGroups = materials->nGroups;
  nR = mesh->rCornerCent.size();
  nZ = mesh->zCornerCent.size();
  nGroupUnknowns = 5*(nZ*nR) + 2*nZ + 2*nR;
  nUnknowns = energyGroups*nGroupUnknowns;

  // initialize size of linear system
  A.resize(nUnknowns,nUnknowns);
  A.reserve(5*nUnknowns+nUnknowns/5);
  x.setZero(nUnknowns);
  xPast.setZero(nUnknowns);
  b.setZero(nUnknowns);

};

//==============================================================================

//==============================================================================
/// Form a portion of the linear system that belongs to SGQD 
/// @param [in] SGQD quasidiffusion energy group to build portion of linear 
///   for
void QDSolver::formLinearSystem(SingleGroupQD * SGQD)	      
{
  int iEq = SGQD->energyGroup*nGroupUnknowns;

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {

      // apply zeroth moment equation
      assertZerothMoment(iR,iZ,iEq,SGQD->energyGroup,SGQD);
      iEq = iEq + 1;

      // north face
      if (iZ == 0)
      {
        // if on the boundary, assert boundary conditions
        assertNFluxBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
        assertNCurrentBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      } else
      {
        // otherwise assert first moment balance on north face
        assertFirstMomentNorth(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      }

      // south face
      if (iZ == mesh->dzsCorner.size()-1)
      {
        // if on the boundary, assert boundary conditions
        assertSFluxBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
        assertSCurrentBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      } else
      {
        // otherwise assert first moment balance on south face
        assertFirstMomentSouth(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      }

      // west face
      if (iR == 0)
      {
        // if on the boundary, assert boundary conditions
        assertWFluxBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
        assertWCurrentBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      } else
      {
        // otherwise assert first moment balance on west face
        assertFirstMomentWest(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      }

      // east face
      if (iR == mesh->drsCorner.size()-1)
      {
        // if on the boundary, assert boundary conditions
        assertEFluxBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
        assertECurrentBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      } else
      {
        // otherwise assert first moment balance on north face
        assertFirstMomentEast(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      }
    }
  }
};

//==============================================================================

//==============================================================================
/// Compute the solution, x, to Ax = b.
void QDSolver::solve()
{
 // Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
 // solver.preconditioner().setDroptol(0.001);
 // solver.compute(A);
 // x = solver.solve(b);
  
  Eigen::SparseLU<Eigen::SparseMatrix<double>,\
    Eigen::COLAMDOrdering<int> > solverLU;
  A.makeCompressed();
  solverLU.compute(A);
  x = solverLU.solve(b);
}
//==============================================================================

//==============================================================================
/// Assert the zeroth moment equation for cell (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert equation for
void QDSolver::assertZerothMoment(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  vector<int> indices;
  vector<double> geoParams = calcGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->neutV(energyGroup);
  double sigT = materials->sigT(iZ,iR,energyGroup);
  double groupSourceCoeff;

  // populate entries representing sources from scattering and fission in 
  // this and other energy groups
  for (int iGroup = 0; iGroup < materials->nGroups; ++iGroup)
  {
    indices = getIndices(iR,iZ,iGroup);
    groupSourceCoeff = calcScatterAndFissionCoeff(iZ,iR,energyGroup,iGroup);
    A.insert(iEq,indices[iCF]) = -geoParams[iCF] * groupSourceCoeff;
  }

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ,energyGroup);

  A.coeffRef(iEq,indices[iCF]) += geoParams[iCF] * ((1/(v*deltaT)) + sigT);

  A.insert(iEq,indices[iWC]) = -geoParams[iWF];

  A.insert(iEq,indices[iEC]) = geoParams[iEF];

  A.insert(iEq,indices[iNC]) = -geoParams[iNF];

  A.insert(iEq,indices[iSC]) = geoParams[iSF];

  // formulate RHS entry
  b(iEq) = geoParams[iCF]*\
    ( (xPast(indices[iCF])/(v*deltaT)) + SGQD->q(iZ,iR));
};
//==============================================================================

//==============================================================================
/// Assert the first moment equation on the south face of cell (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert equation for
void QDSolver::assertFirstMomentSouth(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  vector<int> indices;
  vector<double> geoParams = calcGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->neutV(energyGroup);
  double sigT = materials->sigT(iZ,iR,energyGroup);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zUp-zAvg;

  // get local Eddington factors
  ErrL = SGQD->Err(iZ,iR);
  EzzL = SGQD->Ezz(iZ,iR);
  ErzL = SGQD->Erz(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ,energyGroup);
  A.insert(iEq,indices[iSC]) = ((1/(v*deltaT)) + sigT);

  A.insert(iEq,indices[iSF]) = EzzL/deltaZ;

  A.insert(iEq,indices[iCF]) = -EzzL/deltaZ;

  A.insert(iEq,indices[iWF]) = -(rDown*ErzL/(rAvg*deltaR));

  A.insert(iEq,indices[iEF]) = rUp*ErzL/(rAvg*deltaR);
  
  // formulate RHS entry
  b(iEq) =(xPast(indices[iSC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Assert the first moment equation on the north face of cell (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert equation for
void QDSolver::assertFirstMomentNorth(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  vector<int> indices;
  vector<double> geoParams = calcGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->neutV(energyGroup);
  double sigT = materials->sigT(iZ,iR,energyGroup);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL;  

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rDown; deltaZ = zAvg-zDown;

  // get local Eddington factors
  ErrL = SGQD->Err(iZ,iR);
  EzzL = SGQD->Ezz(iZ,iR);
  ErzL = SGQD->Erz(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ,energyGroup);


  A.insert(iEq,indices[iNC]) = ((1/(v*deltaT)) + sigT);

  A.insert(iEq,indices[iNF]) = -EzzL/deltaZ;

  A.insert(iEq,indices[iCF]) = EzzL/deltaZ;

  A.insert(iEq,indices[iWF]) = -(rDown*ErzL/(rAvg*deltaR));

  A.insert(iEq,indices[iEF]) = rUp*ErzL/(rAvg*deltaR);

  // formulate RHS entry
  b(iEq) =(xPast(indices[iNC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Assert the first moment equation on the west face of cell (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert equation for
void QDSolver::assertFirstMomentWest(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  vector<int> indices;
  vector<double> geoParams = calcGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->neutV(energyGroup);
  double sigT = materials->sigT(iZ,iR,energyGroup);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL;  
  double hCent,hDown;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rAvg-rDown; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg,SGQD);
  hDown = calcIntegratingFactor(iR,iZ,rDown,SGQD);

  // get local Eddington factors
  ErrL = SGQD->Err(iZ,iR);
  EzzL = SGQD->Ezz(iZ,iR);
  ErzL = SGQD->Erz(iZ,iR);
  
  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ,energyGroup);
  A.insert(iEq,indices[iWC]) = ((1/(v*deltaT)) + sigT);

  A.insert(iEq,indices[iSF]) = ErzL/deltaZ;

  A.insert(iEq,indices[iNF]) = -ErzL/deltaZ;

  A.insert(iEq,indices[iCF]) = hCent*ErrL/(hDown*deltaR);

  A.insert(iEq,indices[iWF]) = -hDown*ErrL/(hDown*deltaR);

  // formulate RHS entry
  b(iEq) =(xPast(indices[iWC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Assert the first moment equation on the west face of cell (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert equation for
void QDSolver::assertFirstMomentEast(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  vector<int> indices;
  vector<double> geoParams = calcGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->neutV(energyGroup);
  double sigT = materials->sigT(iZ,iR,energyGroup);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL;  
  double hCent,hUp;

  // calculate geometric values
  rUp = mesh->rCornerEdge(iR+1); rDown = mesh->rCornerEdge(iR);
  zUp = mesh->zCornerEdge(iZ+1); zDown = mesh->zCornerEdge(iZ);
  rAvg = calcVolAvgR(rDown,rUp); zAvg = (zUp + zDown)/2;
  deltaR = rUp-rAvg; deltaZ = zUp-zDown;
  hCent = calcIntegratingFactor(iR,iZ,rAvg,SGQD);
  hUp = calcIntegratingFactor(iR,iZ,rUp,SGQD);

  // get local Eddington factors
  ErrL = SGQD->Err(iZ,iR);
  EzzL = SGQD->Ezz(iZ,iR);
  ErzL = SGQD->Erz(iZ,iR);

  // populate entries representing streaming and reaction terms
  indices = getIndices(iR,iZ,energyGroup);
  A.insert(iEq,indices[iEC]) = ((1/(v*deltaT)) + sigT);

  A.insert(iEq,indices[iSF]) = ErzL/deltaZ;

  A.insert(iEq,indices[iNF]) = -ErzL/deltaZ;

  A.insert(iEq,indices[iCF]) = -hCent*ErrL/(hUp*deltaR);

  A.insert(iEq,indices[iEF]) = hUp*ErrL/(hUp*deltaR);
  
  // formulate RHS entry
  b(iEq) =(xPast(indices[iEC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertNFluxBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  int nIndex = NFluxIndex(iR,iZ,energyGroup);
  
  A.insert(iEq,nIndex) = 1.0;
  b(iEq) = SGQD->nFluxBC(iR);
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertSFluxBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  int sIndex = SFluxIndex(iR,iZ,energyGroup);
  
  A.insert(iEq,sIndex) = 1.0;
  b(iEq) = SGQD->sFluxBC(iR);
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the west face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertWFluxBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  int wIndex = WFluxIndex(iR,iZ,energyGroup);
  
  A.insert(iEq,wIndex) = 1.0;
  b(iEq) = SGQD->wFluxBC(iZ);
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertEFluxBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  int eIndex = EFluxIndex(iR,iZ,energyGroup);
  
  A.insert(iEq,eIndex) = 1.0;
  b(iEq) = SGQD->eFluxBC(iZ);
};
//==============================================================================

//==============================================================================
/// Assert the current boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertNCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  int nIndex = NCurrentIndex(iR,iZ,energyGroup);
  
  A.insert(iEq,nIndex) = 1.0;
  b(iEq) = SGQD->nCurrentZBC(iR);
};
//==============================================================================

//==============================================================================
/// Assert the current boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertSCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  int sIndex = SCurrentIndex(iR,iZ,energyGroup);
  
  A.insert(iEq,sIndex) = 1.0;
  b(iEq) = SGQD->sCurrentZBC(iR);
};
//==============================================================================

//==============================================================================
/// Assert the current boundary condition on the west face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertWCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  int wIndex = WCurrentIndex(iR,iZ,energyGroup);
  
  A.insert(iEq,wIndex) = 1.0;
  b(iEq) = SGQD->wCurrentRBC(iZ);
};
//==============================================================================

//==============================================================================
/// Assert the current boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertECurrentBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  int eIndex = ECurrentIndex(iR,iZ,energyGroup);
  
  A.insert(iEq,eIndex) = 1.0;
  b(iEq) = SGQD->eCurrentRBC(iZ);
};
//==============================================================================

//==============================================================================
/// Calculate multigroup source coefficient for cell at (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] energyGroup energy group of sourcing group
double QDSolver::calcScatterAndFissionCoeff(int iR,int iZ,int toEnergyGroup,\
  int fromEnergyGroup)
{

  double localSigF,localNu,localChiP,localSigS,sourceCoefficient;

  localSigF = materials->sigF(iZ,iR,fromEnergyGroup);
  localNu = materials->nu(iZ,iR);
  localChiP = materials->chiP(iZ,iR,toEnergyGroup);
  localSigS = materials->sigS(iZ,iR,fromEnergyGroup,toEnergyGroup);
  sourceCoefficient = localSigS + localChiP * localNu * localSigF;

  return sourceCoefficient;
};
//==============================================================================

//==============================================================================
/// Form a portion of the linear system that belongs to SGQD 
/// @param [in] SGQD quasidiffusion energy group to build portion of linear 
///   for
double QDSolver::calcIntegratingFactor(int iR,int iZ,double rEval,\
  SingleGroupQD * SGQD)	      
{
  double EzzL,ErrL,ErzL,G,rUp,rDown,rAvg,g0,g1,ratio,hEval;
  int p;
  
  // get local Eddington factors 
  ErrL = SGQD->Err(iZ,iR);
  EzzL = SGQD->Ezz(iZ,iR);
  ErzL = SGQD->Erz(iZ,iR);

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

};

//==============================================================================

//==============================================================================
/// Return global index of average flux at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for average flux in cell at (iR,iZ)
int QDSolver::CFluxIndex(int iR,int iZ,int energyGroup)
{
  int offset = energyGroup*nGroupUnknowns;
  int index = offset + iR + nR * (iZ);
  return index;
};

//==============================================================================

//==============================================================================
/// Return global index of west face flux at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for west face flux in cell at (iR,iZ)
int QDSolver::WFluxIndex(int iR,int iZ,int energyGroup)
{
  int offset = energyGroup*nGroupUnknowns + nR*nZ;
  int westIndex = offset + (iR) + (nR + 1) * iZ;

  return westIndex;
};

//==============================================================================

//==============================================================================
/// Return global index of east face flux at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for east face flux in cell at (iR,iZ)
int QDSolver::EFluxIndex(int iR,int iZ,int energyGroup)
{
  int offset = energyGroup*nGroupUnknowns + nR*nZ;
  int eastIndex = offset + (iR + 1) + (nR + 1) * iZ;

  return eastIndex;
};

//==============================================================================

//==============================================================================
/// Return global index of north face flux at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for north face flux in cell at (iR,iZ)
int QDSolver::NFluxIndex(int iR,int iZ,int energyGroup)
{
  int offset = energyGroup*nGroupUnknowns + nR*nZ + (nR+1)*nZ;
  int northIndex = offset + (iZ) + (nZ + 1) * iR;

  return northIndex;
};

//==============================================================================

//==============================================================================
/// Return global index of south face flux at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for south face flux in cell at (iR,iZ)
int QDSolver::SFluxIndex(int iR,int iZ,int energyGroup)
{
  int offset = energyGroup*nGroupUnknowns + nR*nZ + (nR+1)*nZ;
  int southIndex = offset + (iZ + 1) + (nZ + 1) * iR;

  return southIndex;
};

//==============================================================================

//==============================================================================
/// Return global index of west face current at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for west face current in cell at (iR,iZ)
int QDSolver::WCurrentIndex(int iR,int iZ,int energyGroup)
{
  int offset = energyGroup*nGroupUnknowns + nR*nZ + (nR+1)*nZ + (nZ+1)*nR;
  int westIndex = offset + (iR) + (nR + 1) * iZ;

  return westIndex;
};

//==============================================================================

//==============================================================================
/// Return global index of east face current at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for east face current in cell at (iR,iZ)
int QDSolver::ECurrentIndex(int iR,int iZ,int energyGroup)
{
  int offset = energyGroup*nGroupUnknowns + nR*nZ + (nR+1)*nZ + (nZ+1)*nR;
  int eastIndex = offset + (iR + 1) + (nR + 1) * iZ;

  return eastIndex;
};

//==============================================================================

//==============================================================================
/// Return global index of north face current at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for north face current in cell at (iR,iZ)
int QDSolver::NCurrentIndex(int iR,int iZ,int energyGroup)
{
  int offset = energyGroup*nGroupUnknowns + nR*nZ + 2*(nR+1)*nZ + (nZ+1)*nR;
  int northIndex = offset + (iZ) + (nZ + 1) * iR;

  return northIndex;
};

//==============================================================================

//==============================================================================
/// Return global index of south face current at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for south face current in cell at (iR,iZ)
int QDSolver::SCurrentIndex(int iR,int iZ,int energyGroup)
{
  int offset = energyGroup*nGroupUnknowns + nR*nZ + 2*(nR+1)*nZ + (nZ+1)*nR;
  int southIndex = offset + (iZ + 1) + (nZ + 1) * iR;

  return southIndex;
};

//==============================================================================

//==============================================================================
/// Return global index of south face current at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for south face current in cell at (iR,iZ)
vector<int> QDSolver::getIndices(int iR,int iZ,int energyGroup)
{
  vector<int> indices;
  indices.push_back(CFluxIndex(iR,iZ,energyGroup));
  indices.push_back(WFluxIndex(iR,iZ,energyGroup));
  indices.push_back(EFluxIndex(iR,iZ,energyGroup));
  indices.push_back(NFluxIndex(iR,iZ,energyGroup));
  indices.push_back(SFluxIndex(iR,iZ,energyGroup));
  indices.push_back(WCurrentIndex(iR,iZ,energyGroup));
  indices.push_back(ECurrentIndex(iR,iZ,energyGroup));
  indices.push_back(NCurrentIndex(iR,iZ,energyGroup));
  indices.push_back(SCurrentIndex(iR,iZ,energyGroup));

  return indices;
};

//==============================================================================


//==============================================================================
/// Return geometry parameters for the cell located at (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] gParams vector containing volume and surfaces areas of the 
///   west, east, north, and south faces, in that order.
vector<double> QDSolver::calcGeoParams(int iR,int iZ)
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
double QDSolver::calcVolAvgR(double rDown,double rUp)
{
  // calculate volume-averaged radius
  double volAvgR = (2.0/3.0)*(pow(rUp,3) - pow(rDown,3))/(pow(rUp,2)-pow(rDown,2));

  return volAvgR;
};
//==============================================================================

//==============================================================================
/// Extract cell average values from solution vector and store
/// @param [in] SGQD single group quasidiffusion object to get flux for
void QDSolver::getFlux(SingleGroupQD * SGQD)
{
  vector<int> indices;

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      indices = getIndices(iR,iZ,SGQD->energyGroup);
      SGQD->sFlux(iZ,iR) = x(indices[iCF]);
    }
  }  

};
//==============================================================================

//==============================================================================
/// Extract cell average values from solution vector and store
/// @param [in] SGQD single group quasidiffusion object to get flux for
Eigen::VectorXd QDSolver::getSolutionVector(SingleGroupQD * SGQD)
{
  Eigen::VectorXd solVector(energyGroups*nGroupUnknowns);
  solVector.setZero();
  vector<int> indices;

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      indices = getIndices(iR,iZ,SGQD->energyGroup);
      solVector(indices[iCF]) = SGQD->sFlux(iZ,iR);

      solVector(indices[iWF]) = SGQD->sFluxR(iZ,iR);
      solVector(indices[iEF]) = SGQD->sFluxR(iZ,iR+1);
      solVector(indices[iNF]) = SGQD->sFluxZ(iZ,iR);
      solVector(indices[iSF]) = SGQD->sFluxZ(iZ+1,iR);

      solVector(indices[iWC]) = SGQD->currentR(iZ,iR);
      solVector(indices[iEC]) = SGQD->currentR(iZ,iR+1);
      solVector(indices[iNC]) = SGQD->currentZ(iZ,iR);
      solVector(indices[iSC]) = SGQD->currentZ(iZ+1,iR);
    }
  }  

  return solVector;
};
//==============================================================================








