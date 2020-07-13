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
  
  // Check for optional input parameters 
  checkOptionalParameters();

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
bool MultilevelCoupling::solveOneStep()
{

  Eigen::VectorXd xCurrentIter, xLastIter, residualVector;
  vector<double> residual;
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
    cout << "MGT->MGQD->MPQD Residual: " << residual[0] << endl;
    cout << endl;
    if (residual[0] < mpqd->epsMPQD\
        and residual[1] < mpqd->epsMPQD\
        and eddingtonConverged) 
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
bool MultilevelCoupling::initialSolve()
{

  Eigen::VectorXd xCurrentIter, xLastMGLOQDIter,xLastELOTIter,residualVector;
  vector<double> lastResidualELOT, residualMGLOQD = {1,1,1},\
    residualELOT = {1,1,1};
  bool eddingtonConverged;
  bool convergedMGLOQD=false, convergedELOT=false;
  int itersMGLOQD = 0, itersELOT = 0;

  while (not convergedMGLOQD){

    /////////////////////
    // MGLOQD SOLUTION //
    /////////////////////

    // Solve MGLOQD problem
    solveMGLOQD();

    // Get group fluxes to use in group collapse
    mgqd->getFluxes();

    // Store last iterate of ELOT solution used in MGLOQD level
    xLastMGLOQDIter = mpqd->x;
  
    cout << endl;
    
    while (not convergedELOT) {

      ///////////////////
      // ELOT SOLUTION //
      ///////////////////

      // Calculate collapsed nuclear data
      MGQDToMPQD->collapseNuclearData();

      // Store last iterate of ELOT solution used in ELOT level
      xLastELOTIter = mpqd->x;

      // Solve ELOT problem
      solveELOT();

      // Store newest iterate 
      xCurrentIter = mpqd->x;
     
      // Store last residual   
      lastResidualELOT = residualELOT; 

      // Repeat until ||xCurrentIter - xLastIter|| < eps 
      residualELOT = MGQDToMPQD->calcResidual(xLastELOTIter,xCurrentIter);
      cout << "        ";
      cout << "ELOT Residual: " << residualELOT[0]; 
      cout << ", " << residualELOT[1] << endl;
      itersELOT++;
      
      // Calculate collapnuclear data at new temperature
      mats->updateTemperature(mpqd->heat->returnCurrentTemp());

      // Check if residuals are too big 
     // if (residualELOT[0] > maxResidual or\
     //     residualELOT[1] > maxResidual) 
     // {
      if (residualELOT[0] > lastResidualELOT[0] or\
          residualELOT[1] > lastResidualELOT[1]) 
      {
        // Jump back to MGLOQD level
        break;
      }

      // Check converge criteria 
      if (eps(residualMGLOQD[0], relaxTolELOT) > residualELOT[0] and\
          eps(residualMGLOQD[1], relaxTolELOT) > residualELOT[1]) 
      {
        // ELOT solve is converged
        convergedELOT = true;
      }

    } // ELOT
      
    cout << endl;

    // Reset ELOT residual
    convergedELOT = false; 

    // Store one group fluxes and DNPs for MGHOT and MGLOQD sources 
    mpqd->ggqd->GGSolver->getFlux();
    mpqd->mgdnp->getCumulativeDNPDecaySource();

    // Calculate MGLOQD residual 
    residualMGLOQD = MGQDToMPQD->calcResidual(xLastMGLOQDIter,xCurrentIter);
    cout << "    ";
    cout << "MGLOQD Residual: " << residualMGLOQD[0];
    cout << ", " << residualMGLOQD[1] << endl;
    itersMGLOQD++;

    // Check converge criteria 
    if (eps(mpqd->epsMPQD) > residualMGLOQD[0] and\
        eps(mpqd->epsMPQD) > residualMGLOQD[1])
    { 
      convergedMGLOQD = true;
    }

  } // MGLOQD
    
  cout << endl;
  cout << "MGLOQD iterations: " << itersMGLOQD << endl;
  cout << "ELOT iterations: " << itersELOT << endl;

  return true;

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
bool MultilevelCoupling::solveOneStepLagged()
{

  Eigen::VectorXd xCurrentIter, xLastMGHOTIter, xLastMGLOQDIter,xLastELOTIter,\
    residualVector;
  vector<double> lastResidualELOT, lastResidualMGLOQD, residualMGHOT = {1,1,1},\
    residualMGLOQD = {1,1,1}, residualELOT = {1,1,1};
  bool eddingtonConverged;
  bool convergedMGHOT=false, convergedMGLOQD=false, convergedELOT=false;
  int itersMGHOT = 0, itersMGLOQD = 0, itersELOT = 0;

  while (not convergedMGHOT){ 

    ////////////////////
    // MGHOT SOLUTION // 
    ////////////////////

    // Solve MGHOT problem
    solveMGHOT();

    // Calculate Eddington factors for MGQD problem
    eddingtonConverged = MGTToMGQD->calcEddingtonFactors();

    // Calculate BCs for MGQD problem 
    MGTToMGQD->calcBCs();

    // Store last iterate of ELOT solution used in MGHOT level
    xLastMGHOTIter = mpqd->x;
    
    while (not convergedMGLOQD){

      /////////////////////
      // MGLOQD SOLUTION //
      /////////////////////

      // Solve MGLOQD problem
      solveMGLOQD();

      // Get group fluxes to use in group collapse
      mgqd->getFluxes();
      
      // Store last iterate of ELOT solution used in MGLOQD level
      xLastMGLOQDIter = mpqd->x;
    
      while (not convergedELOT) {
        
        ///////////////////
        // ELOT SOLUTION //
        ///////////////////

        // Calculate collapsed nuclear data
        MGQDToMPQD->collapseNuclearData();
        
        // Store last iterate of ELOT solution used in ELOT level
        xLastELOTIter = mpqd->x;
      
        // Solve ELOT problem
        solveELOT();

        // Store newest iterate 
        xCurrentIter = mpqd->x;

        lastResidualELOT = residualELOT; 

        // Calculate and print ELOT residual  
        residualELOT = MGQDToMPQD->calcResidual(xLastELOTIter,xCurrentIter);
        cout << "        ";
        cout << "ELOT Residual: " << residualELOT[0]; 
        cout << ", " << residualELOT[1] << endl;
        itersELOT++;

        // Calculate collapsed nuclear data at new temperature
        mats->updateTemperature(mpqd->heat->returnCurrentTemp());

        // Check if residuals are too big or if the residuals have increased
        // from the last MGLOQD residual 
//        if (residualELOT[0] > maxResidual or\
//            residualELOT[1] > maxResidual or\
//            residualELOT[0] > lastResidualELOT[0] or\
//            residualELOT[1] > lastResidualELOT[1]) 
//        {
          if (residualELOT[0]/lastResidualELOT[0] > maxResidual or\
              residualELOT[1]/lastResidualELOT[1] > maxResidual) 
          {
          // Jump back to MGLOQD level
          break;
        }
       
        // Check converge criteria 
        if (eps(residualMGLOQD[0], relaxTolELOT) > residualELOT[0] and\
            eps(residualMGLOQD[1], relaxTolELOT) > residualELOT[1]) 
        {
          convergedELOT = true;
        }

      } // ELOT
    
      // Reset convergence indicator
      convergedELOT = false; 
      
      // Store one group fluxes and DNPs for MGHOT and MGLOQD sources 
      mpqd->ggqd->GGSolver->getFlux();
      mpqd->mgdnp->getCumulativeDNPDecaySource();
    
      lastResidualMGLOQD = residualMGLOQD;
 
      // Calculate and print MGLOQD residual 
      residualMGLOQD = MGQDToMPQD->calcResidual(xLastMGLOQDIter,xCurrentIter);
      cout << endl;
      cout << "    ";
      cout << "MGLOQD Residual: " << residualMGLOQD[0];
      cout << ", " << residualMGLOQD[1] << endl;
      cout << endl;
      itersMGLOQD++;

      // Check if residuals are too big or if the residuals have increased
      // from the last MGLOQD residual 
     // if (residualELOT[0] > maxResidual or\
     //     residualELOT[1] > maxResidual or\
     //     residualELOT[0] > lastResidualELOT[0] or\
     //     residualELOT[1] > lastResidualELOT[1]) 
     // {
      if (residualELOT[0]/lastResidualELOT[0] > maxResidual or\
          residualELOT[1]/lastResidualELOT[1] > maxResidual) 
      {
        // Jump back to MGHOT level
        break;
      }


      // Check converge criteria 
      //if (epsELOT(mpqd->epsMPQD) > residualMGLOQD[0] and\
          epsELOT(mpqd->epsMPQD) > residualMGLOQD[1])
      if (eps(residualMGHOT[0],relaxTolMGLOQD) > residualMGLOQD[0] and\
          eps(residualMGHOT[1],relaxTolMGLOQD) > residualMGLOQD[1])
      { 
        convergedMGLOQD = true;
      }

    } // MGLOQD
    
    // Reset convergence indicator
    convergedMGLOQD = false;
     
    // Calculate and print MGHOT residual 
    residualMGHOT = MGQDToMPQD->calcResidual(xLastMGHOTIter,xCurrentIter);
    cout << "MGHOT Residual: " << residualMGHOT[0];
    cout << ", " << residualMGHOT[1] << endl;
    cout << endl;
    itersMGHOT++;
      
    // Check converge criteria 
    if (eps(mpqd->epsMPQD) > residualMGHOT[0] and\
        eps(mpqd->epsMPQD) > residualMGHOT[1]) 
    {
      convergedMGHOT = true;
    }

  } //MGHOT
    
  cout << endl;
    
  cout << "MGHOT iterations: " << itersMGHOT << endl;
  cout << "MGLOQD iterations: " << itersMGLOQD << endl;
  cout << "ELOT iterations: " << itersELOT << endl;

  return true;

};
//==============================================================================

//==============================================================================
/// Dynamically determines a convergence threshold 
///
double MultilevelCoupling::eps(double residual, double relaxationTolerance)
{

  double eps;
  
  if (residual*relaxationTolerance > mpqd->epsMPQD)
    eps = residual*relaxationTolerance;
  else
    eps = mpqd->epsMPQD; 
  
  return eps;

};
//==============================================================================


//==============================================================================
/// Perform a solve at the MGHOT level 
///
void MultilevelCoupling::solveMGHOT()
{

    // Calculate transport sources
    mgt->calcSources();

    // Calculate transport alphas
    mgt->calcAlphas();
  
    // Solve starting angle transport problem
    mgt->solveStartAngles();

    // Solve all angle transport problem
    mgt->solveSCBs();

};
//==============================================================================

//==============================================================================
/// Perform a solve at the MGLOQD level 
///
void MultilevelCoupling::solveMGLOQD()
{

  // Build flux system
  mgqd->buildLinearSystem();

  // Solve flux system
  mgqd->solveLinearSystem();

  // Build neutron current system
  mgqd->buildBackCalcSystem();

  // Solve neutron current system
  mgqd->backCalculateCurrent();

};
//==============================================================================

//==============================================================================
/// Perform a solve at the ELOT level 
///
void MultilevelCoupling::solveELOT()
{

  // Build ELOT system
  mpqd->buildLinearSystem();

  // Solve ELOT system
  mpqd->solveLinearSystem();

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
  //MGQDToMPQD->solveOneStep();
  initialSolve();
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

//==============================================================================
/// Read in optional parameters that might be specified in the input 
///
void MultilevelCoupling::checkOptionalParameters()
{

  // Check for relaxTolELOT specification..
  if ((*input)["parameters"]["relaxTolELOT"])
    relaxTolELOT=(*input)["parameters"]["relaxTolELOT"].as<double>();

  // Check for relaxTolMGLOQD specification.
  if ((*input)["parameters"]["relaxTolELOT"])
    relaxTolMGLOQD=(*input)["parameters"]["relaxTolMGLOQD"].as<double>();

  // Check for maxResidual specification.
  if ((*input)["parameters"]["maxResidual"])
    maxResidual=(*input)["parameters"]["maxResidual"].as<double>();

};
//==============================================================================
