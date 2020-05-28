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
  bool converged;

  // Initialize solve 
  cout << "About to enter MGQDToMPQD::solveOneStep" << endl;
  converged = MGQDToMPQD->solveOneStep();

  cout << "Initial solve completed." << endl;
  cout << endl;

  for (int iStep = 0; iStep < 100; iStep++)
  {
    // Calculate transport sources
    mgt->calcSources();

    cout << "Calculated transport sources." << endl;

    // Calculate transport alphas
    mgt->calcAlphas();
  
    cout << "Calculated alphas." << endl;
    cout << endl;

    // Solve starting angle transport problem
    mgt->solveStartAngles();

    cout << "Solved for starting angles." << endl;
    cout << endl;

    // Solve all angle transport problem
    mgt->solveSCBs();

    cout << "Solved for all angles." << endl;
    cout << endl;
    
    // Calculate Eddington factors for MGQD problem
    MGTToMGQD->calcEddingtonFactors();

    cout << "Calculated Eddington factors." << endl;
    cout << endl;
    
    // Calculate BCs for MGQD problem 
    MGTToMGQD->calcBCs();

    cout << "Calculated MGQD boundary conditions." << endl;
    cout << endl;
    
    // Store xLastIter = MPQD->x
    xLastIter = mpqd->x;

    cout << "Stored xLastIter." << endl;
    cout << endl;
    
    // Solve MGQD and MPQD problem to convergence
    MGQDToMPQD->solveOneStep();

    cout << "Solved coupled MGQD and MPQD problem." << endl;
    cout << endl;
    
    // Store xCurrentIter = MPQD->x
    xCurrentIter = mpqd->x;

    cout << "Stored xCurrentIter." << endl;
    cout << endl;
    
    // Repeat until ||xCurrentIter - xLastIter|| < eps 
    residualVector = xCurrentIter - xLastIter;
    residual = (1.0/xCurrentIter.size())*residualVector.norm();
    cout << "Multilevel Residual: " << residual << endl;
    if (residual < 1E-10) 
    {
      cout << "solve converged!" << endl;
      mgt->writeFluxes();
      break;
    }
  }


};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
void MultilevelCoupling::solveTransient()
{

};
//==============================================================================

