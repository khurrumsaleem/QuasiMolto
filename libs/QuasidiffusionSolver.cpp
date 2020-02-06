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

};

//==============================================================================

//==============================================================================
/// Form a portion of the linear system that belongs to SGQD 
///
/// @param [in] SGQD quasidiffusion energy group to build portion of linear 
///   for
void QDSolver::formLinearSystem(SingleGroupQD * SGQD)	      
{

  cout << "linear system formed" << endl;

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
  int index = iR + nR * (iZ);
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
vector<int> QDSolver::indices(int iR,int iZ,int energyGroup)
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





