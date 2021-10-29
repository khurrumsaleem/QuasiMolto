// File: GreyGroupQD.cpp     
// Purpose: Contains information and builds linear system for grey group qd
// Date: April 9, 2020

#include "GreyGroupQD.h"
#include "GreyGroupSolverTransient.h"
#include "GreyGroupSolverSteadyState.h"
#include "MultiPhysicsCoupledQD.h"

using namespace std;

//==============================================================================
/// GreyGroupQD class object constructor
///
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myInput YAML input object for the simulation
/// @param [in] myMPQD multiphysics coupled QD object for the simulation
GreyGroupQD::GreyGroupQD(Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput,\
    MultiPhysicsCoupledQD * myMPQD)
{
  // Variables for reading inputs
  vector<double> inpSFluxPrev0;
  double fluxSum;

  // Assign pointers for materials, mesh, and input objects
  mpqd = myMPQD;
  materials = myMaterials;
  mesh = myMesh;
  input = myInput;

  GGSolver = std::make_shared<GreyGroupSolver>(this,mesh,materials,input);
  GGSolverBase = std::make_shared<GreyGroupSolverSteadyState>(this,mesh,materials,input);

  // initialize Eddington factors to diffusion physics
  double diagValue = 1.0/3.0, offDiagValue = 0.0;

  Err.setConstant(mesh->nZ,mesh->nR,diagValue);
  Ezz.setConstant(mesh->nZ,mesh->nR,diagValue);
  Erz.setConstant(mesh->nZ,mesh->nR,offDiagValue);
  ErrAxial.setConstant(mesh->nZ+1,mesh->nR,diagValue);
  EzzAxial.setConstant(mesh->nZ+1,mesh->nR,diagValue);
  ErzAxial.setConstant(mesh->nZ+1,mesh->nR,offDiagValue);
  ErrRadial.setConstant(mesh->nZ,mesh->nR+1,diagValue);
  EzzRadial.setConstant(mesh->nZ,mesh->nR+1,diagValue);
  ErzRadial.setConstant(mesh->nZ,mesh->nR+1,offDiagValue);
  GL.setConstant(mesh->nZ,mesh->nR,offDiagValue);
  GR.setConstant(mesh->nZ,mesh->nR,offDiagValue);
  g0.setConstant(mesh->nZ,offDiagValue);
  g1.setConstant(mesh->nZ,offDiagValue);

  // initialize previous Eddington factors
  ErrPrev = Err;
  EzzPrev = Ezz;
  ErzPrev = Erz;

  // initialize source 
  q.setOnes(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  q = 0.0*q;

  // initialize flux and current matrices
  sFlux.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  sFluxR.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size()+1);
  sFluxZ.setZero(mesh->zCornerCent.size()+1,mesh->rCornerCent.size());

  // Check for optional input  
  if ((*input)["parameters"]["initial previous flux"])
  {

    // Read in input
    inpSFluxPrev0=(*input)["parameters"]["initial previous flux"]\
                  .as<vector<double>>();

    // Get sum of group fluxes 
    fluxSum = accumulate(inpSFluxPrev0.begin(),inpSFluxPrev0.end(),0.0);   

    sFlux.setConstant(fluxSum);     
    sFluxR.setConstant(fluxSum);     
    sFluxZ.setConstant(fluxSum);    

  }

  currentR.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size()+1);
  currentZ.setZero(mesh->zCornerCent.size()+1,mesh->rCornerCent.size());

  // initialize boundary conditions
  wFluxBC.setZero(mesh->dzsCorner.size());
  eFluxBC.setZero(mesh->dzsCorner.size());
  nFluxBC.setZero(mesh->drsCorner.size());
  sFluxBC.setZero(mesh->drsCorner.size());

  wCurrentRBC.setZero(mesh->dzsCorner.size());
  eCurrentRBC.setZero(mesh->dzsCorner.size());
  nCurrentZBC.setZero(mesh->drsCorner.size());
  sCurrentZBC.setZero(mesh->drsCorner.size());

  eInwardFluxBC.setZero(mesh->dzsCorner.size());
  nInwardFluxBC.setZero(mesh->drsCorner.size());
  sInwardFluxBC.setZero(mesh->drsCorner.size());

  eInwardCurrentBC.setZero(mesh->dzsCorner.size());
  nInwardCurrentBC.setZero(mesh->drsCorner.size());
  sInwardCurrentBC.setZero(mesh->drsCorner.size());

  eOutwardCurrToFluxRatioBC.setZero(mesh->dzsCorner.size());
  nOutwardCurrToFluxRatioBC.setZero(mesh->drsCorner.size());
  sOutwardCurrToFluxRatioBC.setZero(mesh->drsCorner.size());

  eOutwardCurrToFluxRatioInwardWeightedBC.setZero(mesh->dzsCorner.size());
  nOutwardCurrToFluxRatioInwardWeightedBC.setZero(mesh->drsCorner.size());
  sOutwardCurrToFluxRatioInwardWeightedBC.setZero(mesh->drsCorner.size());

  eAbsCurrentBC.setZero(mesh->dzsCorner.size());
  nAbsCurrentBC.setZero(mesh->drsCorner.size());
  sAbsCurrentBC.setZero(mesh->drsCorner.size());

};
//==============================================================================

//==============================================================================
/// Build linear system for transient QD equations
///
void GreyGroupQD::buildLinearSystem()
{

  GGSolver->formLinearSystem(); // Assuming this is the first set of equations

};
//==============================================================================

//==============================================================================
/// Build linear system for steady state QD equations
///
void GreyGroupQD::buildSteadyStateLinearSystem()
{

  GGSolver->formSteadyStateLinearSystem(); // Assuming this is the first set of equations

};
//==============================================================================

/* PETSc functions */

// Steady state

//==============================================================================
/// Build linear system for steady state QD equations
///
void GreyGroupQD::buildSteadyStateLinearSystem_p()
{

  GGSolver->formSteadyStateLinearSystem_p(); // Assuming this is the first set of equations

};
//==============================================================================

// Transient

//==============================================================================
/// Build linear system for transient QD equations
///
void GreyGroupQD::buildLinearSystem_p()
{

  GGSolver->formLinearSystem_p(); // Assuming this is the first set of equations

};
//==============================================================================

/* Auxiliary functions */

//==============================================================================
/// Print BC parameters 
///
void GreyGroupQD::printBCParams()
{
  cout << "eFluxBC: " << endl;
  cout << eFluxBC << endl;
  cout << endl;

  cout << "nFluxBC: " << endl;
  cout << nFluxBC << endl;
  cout << endl;

  cout << "sFluxBC: " << endl;
  cout << sFluxBC << endl;
  cout << endl;

  cout << "eInwardCurrentBC: " << endl;
  cout << eInwardCurrentBC << endl;
  cout << endl;

  cout << "nInwardCurrentBC: " << endl;
  cout << nInwardCurrentBC << endl;
  cout << endl;

  cout << "sInwardCurrentBC: " << endl;
  cout << sInwardCurrentBC << endl;
  cout << endl;

  cout << "eInwardFluxBC: " << endl;
  cout << eInwardFluxBC << endl;
  cout << endl;

  cout << "nInwardFluxBC: " << endl;
  cout << nInwardFluxBC << endl;
  cout << endl;

  cout << "sInwardFluxBC: " << endl;
  cout << sInwardFluxBC << endl;
  cout << endl;

  cout << "eOutwardCurrToFluxRatioBC: " << endl;
  cout << eOutwardCurrToFluxRatioBC << endl;
  cout << endl;

  cout << "nOutwardCurrToFluxRatioBC: " << endl;
  cout << nOutwardCurrToFluxRatioBC << endl;
  cout << endl;

  cout << "sOutwardCurrToFluxRatioBC: " << endl;
  cout << sOutwardCurrToFluxRatioBC << endl;
  cout << endl;

  cout << "eOutwardCurrToFluxRatioInwardWeightedBC: " << endl;
  cout << eOutwardCurrToFluxRatioInwardWeightedBC << endl;
  cout << endl;

  cout << "nOutwardCurrToFluxRatioInwardWeightedBC: " << endl;
  cout << nOutwardCurrToFluxRatioInwardWeightedBC << endl;
  cout << endl;

  cout << "sOutwardCurrToFluxRatioInwardWeightedBC: " << endl;
  cout << sOutwardCurrToFluxRatioInwardWeightedBC << endl;
  cout << endl;
};
//==============================================================================

//==============================================================================
/// Print BC parameters 
///
void GreyGroupQD::printEddingtons()
{

  cout << "Err: " << endl;
  cout << Err << endl;
  cout << endl;

  cout << "Ezz: " << endl;
  cout << Ezz << endl;
  cout << endl;

  cout << "Erz: " << endl;
  cout << Erz << endl;
  cout << endl;

  cout << "ErrAxial: " << endl;
  cout << ErrAxial << endl;
  cout << endl;

  cout << "EzzAxial: " << endl;
  cout << EzzAxial << endl;
  cout << endl;

  cout << "ErzAxial: " << endl;
  cout << ErzAxial << endl;
  cout << endl;

  cout << "ErrRadial: " << endl;
  cout << ErrRadial << endl;
  cout << endl;

  cout << "EzzRadial: " << endl;
  cout << EzzRadial << endl;
  cout << endl;

  cout << "ErzRadial: " << endl;
  cout << ErzRadial << endl;
  cout << endl;

};
//==============================================================================


//==============================================================================
/// Assign pointer to MultiGroupQuasiDiffusion object 
///
void GreyGroupQD::assignMGQDPointer(MultiGroupQD * myMGQD)
{

  mgqd = myMGQD;
  useMGQDSources = true; 

};
//==============================================================================




