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
    if (residual < mpqd->epsMPQD and eddingtonConverged) 
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
bool MultilevelCoupling::solveOneStepLagged()
{

  Eigen::VectorXd xCurrentIter, xLastMGHOTIter, xLastMGLOQDIter,xLastELOTIter,\
    residualVector;
  double residualMGHOT=1,residualMGLOQD=2,residualELOT=3;
  bool eddingtonConverged,convergedMGHOT,convergedMGLOQD,convergedELOT;
  int itersMGHOT = 0, itersMGLOQD = 0, itersELOT = 0;

  while (mpqd->epsMPQD < residualMGHOT){ 

    //////////////////////// 
    // TRANSPORT SOLUTION // 
    //////////////////////// 

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

    // Store last iterate of ELOT solution used in MGHOT level
    xLastMGHOTIter = mpqd->x;
    
    //while (mpqd->epsMPQD < residualMGLOQD){
    while (mpqd->epsMPQD < residualMGLOQD){

      /////////////////////
      // MGLOQD SOLUTION //
      /////////////////////

      // Build flux system
      mgqd->buildLinearSystem();
      
      // Solve flux system
      mgqd->solveLinearSystem();

      // Build neutron current system
      mgqd->buildBackCalcSystem();

      // Solve neutron current system
      mgqd->backCalculateCurrent();

      // Get group fluxes to use in group collapse
      mgqd->getFluxes();
      
      // Store last iterate of ELOT solution used in MGLOQD level
      xLastMGLOQDIter = mpqd->x;
    
      //while (mpqd->epsMPQD < residualELOT){
      while ((residualMGLOQD*1E-5 < residualELOT\
          and mpqd->epsMPQD < residualELOT)){
        
        ///////////////////
        // ELOT SOLUTION //
        ///////////////////

        // Calculate collapsed nuclear data
        MGQDToMPQD->collapseNuclearData();
        
        // Store last iterate of ELOT solution used in ELOT level
        xLastELOTIter = mpqd->x;

        // Build ELOT system
        mpqd->buildLinearSystem();

        // Solve ELOT system
        mpqd->solveLinearSystem();

        // Store newest iterate 
        xCurrentIter = mpqd->x;

        // Repeat until ||xCurrentIter - xLastIter|| < eps 
        residualELOT = MGQDToMPQD->calcResidual(xLastELOTIter,xCurrentIter);
        cout << endl; 
        cout << "        ELOT Residual: " << residualELOT << endl;
        cout << endl;
        itersELOT++;
        
        // Calculate collapsed nuclear data at new temperature
        mats->updateTemperature(mpqd->heat->returnCurrentTemp());

      } // ELOT
     
      // Reset ELOT residual
      residualELOT = 3; 
      
      // Store one group fluxes and DNPs for MGHOT and MGLOQD sources 
      mpqd->ggqd->GGSolver->getFlux();
      mpqd->mgdnp->getCumulativeDNPDecaySource();
     
      // Calculate MGLOQD residual 
      residualMGLOQD = MGQDToMPQD->calcResidual(xLastMGLOQDIter,xCurrentIter);
      cout << endl; 
      cout << "    MGLOQD Residual: " << residualMGLOQD << endl;
      cout << endl;
      itersMGLOQD++;

    } // MGLOQD
    
    // Reset MGLOQD residual
    residualMGLOQD = 2; 
     
    // Calculate MGHOT residual 
    residualMGHOT = MGQDToMPQD->calcResidual(xLastMGHOTIter,xCurrentIter);
    cout << endl; 
    cout << "MGHOT Residual: " << residualMGHOT << endl;
    cout << endl;
    itersMGHOT++;

  } //MGHOT
    
  cout << "MGHOT iterations: " << itersMGHOT << endl;
  cout << "MGLOQD iterations: " << itersMGLOQD << endl;
  cout << "ELOT iterations: " << itersELOT << endl;

  return true;

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
    if(solveOneStepLagged())
    {
      mgqd->updateVarsAfterConvergence(); 
      mpqd->updateVarsAfterConvergence(); 
      mgqd->writeVars();
      mpqd->writeVars(); 
      mats->oneGroupXS->writeVars();
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

