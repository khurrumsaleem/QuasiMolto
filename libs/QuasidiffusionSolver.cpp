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
  int nUnknowns,nCurrentUnknowns;

  // calculate number of unknowns  
  energyGroups = materials->nGroups;
  nR = mesh->rCornerCent.size();
  nZ = mesh->zCornerCent.size();
  nGroupUnknowns = 3*(nZ*nR) + nZ + nR;
  nGroupCurrentUnknowns = 2*(nZ*nR) + nZ + nR;
  nUnknowns = energyGroups*nGroupUnknowns;
  nCurrentUnknowns = energyGroups*nGroupCurrentUnknowns;

  // initialize size of linear system
  A.resize(nUnknowns,nUnknowns);
  A.reserve(3*nUnknowns+nUnknowns/5);
  C.resize(nCurrentUnknowns,nUnknowns);
  C.reserve(4*nCurrentUnknowns);
  x.setZero(nUnknowns);
  xPast.setZero(nUnknowns);
  currPast.setZero(energyGroups*nGroupCurrentUnknowns);
  b.setZero(nUnknowns);
  d.setZero(nCurrentUnknowns);

  checkOptionalParams();
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
        assertNBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      } 

      // south face
      if (iZ == mesh->dzsCorner.size()-1)
      {
        // if on the boundary, assert boundary conditions
        assertSBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      } else
      {
        // otherwise assert first moment balance on south face
        applyAxialBoundary(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      }

      // west face
      if (iR == 0)
      {
        // if on the boundary, assert boundary conditions
        assertWBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      } 

      // east face
      if (iR == mesh->drsCorner.size()-1)
      {
        // if on the boundary, assert boundary conditions
        assertEBC(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      } else
      {
        // otherwise assert first moment balance on north face
        applyRadialBoundary(iR,iZ,iEq,SGQD->energyGroup,SGQD);
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
/// Compute currents
void QDSolver::backCalculateCurrent()
{
  currPast = d + C*x; 
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

  westCurrent(-geoParams[iWF],iR,iZ,iEq,energyGroup,SGQD);
  
  eastCurrent(geoParams[iEF],iR,iZ,iEq,energyGroup,SGQD);

  northCurrent(-geoParams[iNF],iR,iZ,iEq,energyGroup,SGQD);

  southCurrent(geoParams[iSF],iR,iZ,iEq,energyGroup,SGQD);

  // formulate RHS entry
  b(iEq) = b(iEq) + geoParams[iCF]*\
    ( (xPast(indices[iCF])/(v*deltaT)) + SGQD->q(iZ,iR));
};
//==============================================================================

//==============================================================================
void QDSolver::applyRadialBoundary(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  eastCurrent(1,iR,iZ,iEq,energyGroup,SGQD);
  westCurrent(-1,iR+1,iZ,iEq,energyGroup,SGQD);
}
//==============================================================================

//==============================================================================
void QDSolver::applyAxialBoundary(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  northCurrent(1,iR,iZ+1,iEq,energyGroup,SGQD);
  southCurrent(-1,iR,iZ,iEq,energyGroup,SGQD);
}
//==============================================================================

//==============================================================================
/// Enforce coefficients for current on south face
void QDSolver::southCurrent(double coeff,int iR,int iZ,int iEq,int energyGroup,\
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

  coeff = coeff/((1/(v*deltaT))+sigT); 

  A.coeffRef(iEq,indices[iSF]) -= coeff*EzzL/deltaZ;

  A.coeffRef(iEq,indices[iCF]) += coeff*EzzL/deltaZ;

  A.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzL/(rAvg*deltaR));

  A.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzL/(rAvg*deltaR);
  
  // formulate RHS entry
  b(iEq) = b(iEq) - coeff*(currPast(indices[iSC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for current on north face 
void QDSolver::northCurrent(double coeff,int iR,int iZ,int iEq,int energyGroup,\
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
  
  coeff = coeff/((1/(v*deltaT))+sigT); 

  A.coeffRef(iEq,indices[iNF]) += coeff*EzzL/deltaZ;

  A.coeffRef(iEq,indices[iCF]) -= coeff*EzzL/deltaZ;

  A.coeffRef(iEq,indices[iWF]) += coeff*(rDown*ErzL/(rAvg*deltaR));

  A.coeffRef(iEq,indices[iEF]) -= coeff*rUp*ErzL/(rAvg*deltaR);

  // formulate RHS entry
  b(iEq) = b(iEq) - coeff*(currPast(indices[iNC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for current on west face
void QDSolver::westCurrent(double coeff,int iR,int iZ,int iEq,int energyGroup,\
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
  
  coeff = coeff/((1/(v*deltaT))+sigT); 

  A.coeffRef(iEq,indices[iSF]) -= coeff*ErzL/deltaZ;

  A.coeffRef(iEq,indices[iNF]) += coeff*ErzL/deltaZ;

  A.coeffRef(iEq,indices[iCF]) -= coeff*hCent*ErrL/(hDown*deltaR);

  A.coeffRef(iEq,indices[iWF]) += coeff*hDown*ErrL/(hDown*deltaR);

  // formulate RHS entry
  b(iEq) = b(iEq) - coeff*(currPast(indices[iWC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Enforce coefficients for current on east face
void QDSolver::eastCurrent(double coeff,int iR,int iZ,int iEq,int energyGroup,\
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
  
  coeff = coeff/((1/(v*deltaT))+sigT); 

  A.coeffRef(iEq,indices[iSF]) -= coeff*ErzL/deltaZ;

  A.coeffRef(iEq,indices[iNF]) += coeff*ErzL/deltaZ;

  A.coeffRef(iEq,indices[iCF]) += coeff*hCent*ErrL/(hUp*deltaR);

  A.coeffRef(iEq,indices[iEF]) -= coeff*hUp*ErrL/(hUp*deltaR);
  
  // formulate RHS entry
  b(iEq) = b(iEq) - coeff*(currPast(indices[iEC])/(v*deltaT));
};
//==============================================================================

//==============================================================================
/// Form a portion of the current back calc linear system that belongs to SGQD 
/// @param [in] SGQD quasidiffusion energy group to build portion of linear 
///   for
void QDSolver::formBackCalcSystem(SingleGroupQD * SGQD)	      
{
  int iEq = SGQD->energyGroup*nGroupCurrentUnknowns;

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      
      // south face
      calcSouthCurrent(iR,iZ,iEq,SGQD->energyGroup,SGQD);
      iEq = iEq + 1;

      // east face
      calcEastCurrent(iR,iZ,iEq,SGQD->energyGroup,SGQD);
      iEq = iEq + 1;

      // north face
      if (iZ == 0)
      {
        calcNorthCurrent(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      } 

      // west face
      if (iR == 0)
      {
        // if on the boundary, assert boundary conditions
        calcWestCurrent(iR,iZ,iEq,SGQD->energyGroup,SGQD);
        iEq = iEq + 1;
      } 

    }
  }
};

//==============================================================================


//==============================================================================
/// Enforce coefficients to calculate current on south face
void QDSolver::calcSouthCurrent(int iR,int iZ,int iEq,\
  int energyGroup,SingleGroupQD * SGQD)
{
  vector<int> indices;
  vector<double> geoParams = calcGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->neutV(energyGroup);
  double sigT = materials->sigT(iZ,iR,energyGroup);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL,coeff;  

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
void QDSolver::calcNorthCurrent(int iR,int iZ,int iEq,\
  int energyGroup,SingleGroupQD * SGQD)
{
  vector<int> indices;
  vector<double> geoParams = calcGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->neutV(energyGroup);
  double sigT = materials->sigT(iZ,iR,energyGroup);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL,coeff;  

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
void QDSolver::calcWestCurrent(int iR,int iZ,int iEq,\
  int energyGroup,SingleGroupQD * SGQD)
{
  vector<int> indices;
  vector<double> geoParams = calcGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->neutV(energyGroup);
  double sigT = materials->sigT(iZ,iR,energyGroup);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL,coeff;  
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
void QDSolver::calcEastCurrent(int iR,int iZ,int iEq,\
  int energyGroup,SingleGroupQD * SGQD)
{
  vector<int> indices;
  vector<double> geoParams = calcGeoParams(iR,iZ);
  double deltaT = mesh->dt;
  double v = materials->neutV(energyGroup);
  double sigT = materials->sigT(iZ,iR,energyGroup);
  double rUp,rDown,zUp,zDown,rAvg,zAvg,deltaR,deltaZ,ErrL,EzzL,ErzL,coeff;  
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
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertNFluxBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  vector<int> indices = getIndices(iR,iZ,energyGroup);
  
  A.insert(iEq,indices[iNF]) = 1.0;
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
  vector<int> indices = getIndices(iR,iZ,energyGroup);
  
  A.insert(iEq,indices[iSF]) = 1.0;
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
  vector<int> indices = getIndices(iR,iZ,energyGroup);
  
  A.insert(iEq,indices[iWF]) = 1.0;
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
  vector<int> indices = getIndices(iR,iZ,energyGroup);
  
  A.insert(iEq,indices[iEF]) = 1.0;
  b(iEq) = SGQD->eFluxBC(iZ);
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertNCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  northCurrent(1,iR,iZ,iEq,energyGroup,SGQD);
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertSCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  southCurrent(1,iR,iZ,iEq,energyGroup,SGQD);
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the west face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertWCurrentBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  westCurrent(1,iR,iZ,iEq,energyGroup,SGQD);
};
//==============================================================================

//==============================================================================
/// Assert the flux boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertECurrentBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  eastCurrent(1,iR,iZ,iEq,energyGroup,SGQD);
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertNGoldinBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  vector<int> indices = getIndices(iR,iZ,energyGroup);
  double ratio = SGQD->nOutwardCurrToFluxRatioBC(iR);
  double inwardCurrent = SGQD->nInwardCurrentBC(iR);
  double inwardFlux = SGQD->nInwardFluxBC(iR);

  northCurrent(1,iR,iZ,iEq,energyGroup,SGQD);
  A.coeffRef(iEq,indices[iNF]) -= ratio;
  b(iEq) = b(iEq) + (inwardCurrent-ratio*inwardFlux);
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertSGoldinBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  vector<int> indices = getIndices(iR,iZ,energyGroup);
  double ratio = SGQD->sOutwardCurrToFluxRatioBC(iR);
  double inwardCurrent = SGQD->sInwardCurrentBC(iR);
  double inwardFlux = SGQD->sInwardFluxBC(iR);

  southCurrent(1,iR,iZ,iEq,energyGroup,SGQD);
  A.coeffRef(iEq,indices[iSF]) -= ratio;
  b(iEq) = b(iEq) + (inwardCurrent-ratio*inwardFlux);
};
//==============================================================================

//==============================================================================
/// Assert Gol'din's boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertEGoldinBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  vector<int> indices = getIndices(iR,iZ,energyGroup);
  double ratio = SGQD->eOutwardCurrToFluxRatioBC(iR);
  double inwardCurrent = SGQD->eInwardCurrentBC(iR);
  double inwardFlux = SGQD->eInwardFluxBC(iR);

  eastCurrent(1,iR,iZ,iEq,energyGroup,SGQD);
  A.coeffRef(iEq,indices[iEF]) -= ratio;
  b(iEq) = b(iEq) + (inwardCurrent-ratio*inwardFlux);
};
//==============================================================================


//==============================================================================
/// Assert boundary condition on the north face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertNBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  if (reflectingBCs)
    assertNCurrentBC(iR,iZ,iEq,energyGroup,SGQD);
  else if (reflectingBCs)
    assertNGoldinBC(iR,iZ,iEq,energyGroup,SGQD);
  else
    assertNFluxBC(iR,iZ,iEq,energyGroup,SGQD);
};
//==============================================================================

//==============================================================================
/// Assert boundary condition on the south face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertSBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  if (reflectingBCs)
    assertSCurrentBC(iR,iZ,iEq,energyGroup,SGQD);
  else if (goldinBCs)
    assertSGoldinBC(iR,iZ,iEq,energyGroup,SGQD);
  else
    assertSFluxBC(iR,iZ,iEq,energyGroup,SGQD);
};
//==============================================================================

//==============================================================================
/// Assert boundary condition on the west face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertWBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  if (reflectingBCs or goldinBCs)
    assertWCurrentBC(iR,iZ,iEq,energyGroup,SGQD);
  else
    assertWFluxBC(iR,iZ,iEq,energyGroup,SGQD);
};
//==============================================================================

//==============================================================================
/// Assert boundary condition on the east face at location (iR,iZ)
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [in] iEq row to place equation in
/// @param [in] energyGroup energy group to assert boundary condition for
void QDSolver::assertEBC(int iR,int iZ,int iEq,int energyGroup,\
  SingleGroupQD * SGQD)
{
  if (reflectingBCs)
    assertECurrentBC(iR,iZ,iEq,energyGroup,SGQD);
  if (goldinBCs)
    assertEGoldinBC(iR,iZ,iEq,energyGroup,SGQD);
  else
    assertEFluxBC(iR,iZ,iEq,energyGroup,SGQD);
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
/// Return global index of south face current at indices iR and iZ 
/// @param [in] iR radial index of cell
/// @param [in] iZ axial index of cell
/// @param [out] index global index for south face current in cell at (iR,iZ)
vector<int> QDSolver::getIndices(int iR,int iZ,int energyGroup)
{
  
  vector<int> indices,oneGroupIndices;
  
  // Set flux and current offsets according to energy group
  int offsetFlux = energyGroup*nGroupUnknowns;
  int offsetCurr = energyGroup*nGroupCurrentUnknowns;

  // Get indices for a single energy group 
  oneGroupIndices = mesh->getQDCellIndices(iR,iZ);

  // Offset by specified energy group
  indices.push_back(oneGroupIndices[iCF] + offsetFlux);
  indices.push_back(oneGroupIndices[iWF] + offsetFlux);
  indices.push_back(oneGroupIndices[iEF] + offsetFlux);
  indices.push_back(oneGroupIndices[iNF] + offsetFlux);
  indices.push_back(oneGroupIndices[iSF] + offsetFlux);
  indices.push_back(oneGroupIndices[iWC] + offsetCurr);
  indices.push_back(oneGroupIndices[iEC] + offsetCurr);
  indices.push_back(oneGroupIndices[iNC] + offsetCurr);
  indices.push_back(oneGroupIndices[iSC] + offsetCurr);

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
/// Extract flux values from SGQD object, store to solution vector, and return 
/// @param [in] SGQD single group quasidiffusion object to get flux for
Eigen::VectorXd QDSolver::getFluxSolutionVector(SingleGroupQD * SGQD)
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

    }
  }  

  return solVector;
};
//==============================================================================

//==============================================================================
/// Extract current values from SGQD object, store to solution vector, and 
/// return 
/// @param [in] SGQD single group quasidiffusion object to get current for
Eigen::VectorXd QDSolver::getCurrentSolutionVector(SingleGroupQD * SGQD)
{
  Eigen::VectorXd solVector(energyGroups*nGroupCurrentUnknowns);
  solVector.setZero();
  vector<int> indices;

  // loop over spatial mesh
  for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
    {
      indices = getIndices(iR,iZ,SGQD->energyGroup);

      solVector(indices[iWC]) = SGQD->currentR(iZ,iR);
      solVector(indices[iEC]) = SGQD->currentR(iZ,iR+1);
      solVector(indices[iNC]) = SGQD->currentZ(iZ,iR);
      solVector(indices[iSC]) = SGQD->currentZ(iZ+1,iR);
    }
  }  

  return solVector;
};
//==============================================================================

void QDSolver::checkOptionalParams()
{
  string boundaryType;

  // check for optional parameters specified in input file

  if ((*input)["parameters"]["solve type"])
  {

    boundaryType=(*input)["parameters"]["solve type"].as<string>();
    if (boundaryType == "TQD") goldinBCs = true;

  }
  else if ((*input)["parameters"]["mgqd-bcs"])
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
