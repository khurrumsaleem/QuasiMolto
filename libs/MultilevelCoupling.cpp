// File: MultilevelCoupling.cpp     
// Purpose: Solve a coupled MGT, MGQD, MPQD problem. 
// Date: May 27, 2020

#include "MultilevelCoupling.h"

using namespace std;

//==============================================================================
/// MultilevelCoupling class object constructor
///
/// @param [in] myMesh mesh object 
/// @param [in] myMats materials object 
/// @param [in] myInput input object 
/// @param [in] myMPQD MultiPhysicsCoupledQD object 
/// @param [in] myMGPQD multigroup quasidiffusion object 
MultilevelCoupling::MultilevelCoupling(Mesh * myMesh,\
        Materials * myMats,\
        YAML::Node * myInput,\
        MultiGroupTransport * myMGT,\
        MultiGroupQD * myMGQD,\
        MultiPhysicsCoupledQD * myMPQD)

{

  // Assign inputs to their member variables
  mesh = myMesh; 
  mats = myMats;
  input = myInput;
  mgt = myMGT;
  mpqd = myMPQD;
  mgqd = myMGQD;

  // Assign MPQD pointer in multigroup QD solver
  mgt->assignMultiPhysicsCoupledQDPointers(mgqd,mpqd);

  // Assign MPQD pointer in multigroup QD solver
  mgqd->assignMultiPhysicsCoupledQDPointer(mpqd);

  // Assign MGQD pointer in grey group solver
  mpqd->ggqd->assignMGQDPointer(mgqd);

  // Create MGT to MGQD coupling object
  MGTToMGQD = new TransportToQDCoupling(mats,mesh,input,mgt,mgqd);
 
  // Create MGQD to MPQD coupling object
  MGQDToMPQD = new MGQDToMPQDCoupling(mesh,mats,input,mpqd,mgqd);
  

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
bool MultilevelCoupling::solveOneStep()
{

  Eigen::VectorXd xCurrentIter, xLastIter, residualVector;
  double residual;
  bool eddingtonConverged;

  for (int iStep = 0; iStep < 100; iStep++)
  {
    // Calculate transport sources
    mgt->calcSources();

    // Calculate transport alphas
    mgt->calcAlphas();
  
    // Solve starting angle transport problem
    mgt->solveStartAngles();

    // Solve all angle transport problem
    mgt->solveSCBs();

    // Calculate Eddington factors for MGQD problem
    eddingtonConverged = MGTToMGQD->calcEddingtonFactors();

    // Calculate BCs for MGQD problem 
    MGTToMGQD->calcBCs();

    // Store xLastIter = MPQD->x
    xLastIter = mpqd->x;

    // Solve MGQD and MPQD problem to convergence
    MGQDToMPQD->solveOneStep();

    // Store xCurrentIter = MPQD->x
    xCurrentIter = mpqd->x;

    // Repeat until ||xCurrentIter - xLastIter|| < eps 
    residual = MGQDToMPQD->calcResidual(xLastIter,xCurrentIter);
    cout << endl; 
    cout << "MGT->MGQD->MPQD Residual: " << residual << endl;
    cout << endl;
    if (residual < 1E-10 and eddingtonConverged) 
    {
      cout << "Solve converged." << endl;
      mgt->writeFluxes();
      return true;
    }
  }

  return false;

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
void MultilevelCoupling::solveTransient()
{
  // Write mesh info
  mesh->writeVars();

  // Initialize solve 
  cout << "Computing initial solve..." << endl;
  cout << endl;
  MGQDToMPQD->solveOneStep();
  cout << "Initial solve completed." << endl;
  cout << endl;
  
  for (int iTime = 0; iTime < mesh->dts.size(); iTime++)
  {
  cout << "Solve for t = "<< mesh->ts[iTime+1] << endl;
  cout << endl;
    if(solveOneStep())
    {
      mgqd->updateVarsAfterConvergence(); 
      mpqd->updateVarsAfterConvergence(); 
      mgqd->writeVars();
      mpqd->writeVars();
      mesh->advanceOneTimeStep();
    }  
    else 
    {
    cout << "Solve not converged." << endl;
    cout << "Transient aborted." << endl;
    break;
    }
  } 

};
//==============================================================================

