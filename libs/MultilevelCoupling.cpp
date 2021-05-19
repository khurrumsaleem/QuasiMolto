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
      solveELOT(mpqd->x);

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
      if (residualELOT[0]/lastResidualELOT[0] > resetThreshold or\
          residualELOT[1]/lastResidualELOT[0] > resetThreshold) 
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
bool MultilevelCoupling::solveOneStepResidualBalance(bool outputVars)
{

  Eigen::VectorXd xCurrentIter, xLastMGHOTIter, xLastMGLOQDIter,xLastELOTIter,\
    residualVector;
  vector<double> lastResidualELOT, lastResidualMGLOQD, residualMGHOT = {1,1,1},\
    residualMGLOQD = {1,1,1}, residualELOT = {1,1,1};
  bool eddingtonConverged;
  bool convergedMGHOT=false, convergedMGLOQD=false, convergedELOT=false;
  int itersMGHOT = 0, itersMGLOQD = 0, itersELOT = 0;
  vector<int> iters;
  vector<double> tempResMGHOT,tempResMGLOQD,tempResELOT,tempResiduals;
  vector<double> fluxResMGHOT,fluxResMGLOQD,fluxResELOT,fluxResiduals;
 
  // Timing variables 
  double duration,totalDuration = 0.0,elotDuration = 0,\
    mgloqdDuration = 0, mghotDuration = 0;
  clock_t startTime;
  auto begin = chrono::high_resolution_clock::now();

  while (not convergedMGHOT){ 

    ////////////////////
    // MGHOT SOLUTION // 
    ////////////////////

    // Only solve the MGHOT after we've got an estimate for the ELOT solution
    if (itersMGHOT != 0 and not p1Approx)
    {
      // Solve MGHOT problem
      cout << "MGHOT solve...";
      //startTime = clock(); 
      auto begin = chrono::high_resolution_clock::now();
      solveMGHOT();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      //duration = (clock() - startTime)/(double)CLOCKS_PER_SEC;
      mghotDuration = mghotDuration + duration;
      cout << " done. ("<< duration << " seconds)" << endl;
      iters.push_back(3);

      // Calculate Eddington factors for MGQD problem
      cout << "Calculating MGLOQD Eddington factors...";
      eddingtonConverged = MGTToMGQD->calcEddingtonFactors();
      cout << " done."<< endl;

      // Calculate BCs for MGQD problem 
      cout << "Calculating MGLOQD boundary conditions...";
      MGTToMGQD->calcBCs();
      cout << " done."<< endl;
    }

    // Store last iterate of ELOT solution used in MGHOT level
    xLastMGHOTIter = mpqd->x;
    
    while (not convergedMGLOQD){

      /////////////////////
      // MGLOQD SOLUTION //
      /////////////////////

      // Solve MGLOQD problem
      cout << "    ";
      cout << "MGLOQD solve..." << endl;
      //startTime = clock(); 
      auto begin = chrono::high_resolution_clock::now();
      solveMGLOQD();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      //duration = (clock() - startTime)/(double)CLOCKS_PER_SEC;
      mgloqdDuration = mgloqdDuration + duration;
      cout << "    ";
      cout << "MGLOQD solve done. ("<< duration << " seconds)" << endl;
      iters.push_back(2);

      // Get group fluxes to use in group collapse
      mgqd->getFluxes();
      
      // Store last iterate of ELOT solution used in MGLOQD level
      xLastMGLOQDIter = mpqd->x;
    
      while (not convergedELOT) {
        
        ///////////////////
        // ELOT SOLUTION //
        ///////////////////

        // Calculate collapsed nuclear data
        cout << "        ";
        cout << "Collapsing MGLOQD data for ELOT solve...";
        MGQDToMPQD->collapseNuclearData();
        cout << " done."<< endl;
        
        // Store last iterate of ELOT solution used in ELOT level
        cout << "        ";
        cout << "Storing last solution...";
        xLastELOTIter = mpqd->x;
        cout << " done."<< endl;
      
        // Solve ELOT problem
        cout << "        ";
        cout << "ELOT solve..." << endl;
        //startTime = clock(); 
        auto begin = chrono::high_resolution_clock::now();
        solveELOT(mpqd->x);
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        duration = elapsed.count()*1e-9;
        //duration = (clock() - startTime)/(double)CLOCKS_PER_SEC;
        elotDuration = elotDuration + duration;
        cout << "        ";
        cout << "ELOT solve done. ("<< duration << " seconds)" << endl;
        iters.push_back(1);

        // Store newest iterate 
        xCurrentIter = mpqd->x;

        lastResidualELOT = residualELOT; 

        // Calculate and print ELOT residual  
        residualELOT = MGQDToMPQD->calcResidual(xLastELOTIter,xCurrentIter);
        residualELOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastELOTIter);
        cout << "        ";
        cout << "ELOT Residual: " << residualELOT[0]; 
        cout << ", " << residualELOT[1] << endl;

        // Update iterate counters, store residuals
        itersELOT++;
        fluxResELOT.push_back(residualELOT[0]);
        tempResELOT.push_back(residualELOT[1]);

        // Calculate collapsed nuclear data at new temperature
        mats->updateTemperature(mpqd->heat->returnCurrentTemp());

        // Check if residuals are too big or if the residuals have increased
        // from the last MGLOQD residual 
        if (residualELOT[0]/lastResidualELOT[0] > resetThreshold and\
            residualELOT[1]/lastResidualELOT[1] > resetThreshold) 
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
      residualMGLOQD = MGQDToMPQD->calcResidual(xCurrentIter,xLastMGLOQDIter);
      cout << endl;
      cout << "    ";
      cout << "MGLOQD Residual: " << residualMGLOQD[0];
      cout << ", " << residualMGLOQD[1] << endl;
      cout << endl;
        
      // Update iterate counters, store residuals
      itersMGLOQD++;

      fluxResMGLOQD.push_back(residualMGLOQD[0]);
      fluxResMGLOQD.insert(fluxResMGLOQD.end(),\
          fluxResELOT.begin(),fluxResELOT.end());
      fluxResELOT.clear();

      tempResMGLOQD.push_back(residualMGLOQD[1]);
      tempResMGLOQD.insert(tempResMGLOQD.end(),\
          tempResELOT.begin(),tempResELOT.end());
      tempResELOT.clear();

      // Check if residuals are too big or if the residuals have increased
      // from the last MGLOQD residual 
      if (residualMGLOQD[0]/lastResidualMGLOQD[0] > resetThreshold or\
          residualMGLOQD[1]/lastResidualMGLOQD[1] > resetThreshold) 
      {
        // Jump back to MGHOT level
        break;
      }

      // Check converge criteria 
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
    residualMGHOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastMGHOTIter);
    cout << "MGHOT Residual: " << residualMGHOT[0];
    cout << ", " << residualMGHOT[1] << endl;
    cout << endl;
      
    // Update iterate counters, store residuals
    itersMGHOT++;
    
    if (itersMGHOT != 1) fluxResMGHOT.push_back(residualMGHOT[0]);
    fluxResMGHOT.insert(fluxResMGHOT.end(),\
        fluxResMGLOQD.begin(),fluxResMGLOQD.end());
    fluxResMGLOQD.clear();
    
    if (itersMGHOT != 1) tempResMGHOT.push_back(residualMGHOT[1]);
    tempResMGHOT.insert(tempResMGHOT.end(),\
        tempResMGLOQD.begin(),tempResMGLOQD.end());
    tempResMGLOQD.clear();

    // Update persistent iteration record
    fluxResiduals.insert(fluxResiduals.end(),fluxResMGHOT.begin(),\
        fluxResMGHOT.end());
    fluxResMGHOT.clear();
    
    tempResiduals.insert(tempResiduals.end(),tempResMGHOT.begin(),\
        tempResMGHOT.end());
    tempResMGHOT.clear();
      
    // Check converge criteria 
    if (eps(mpqd->epsMPQD) > residualMGHOT[0] and\
        eps(mpqd->epsMPQD) > residualMGHOT[1]) 
    {
      convergedMGHOT = true;
    }

  } //MGHOT
    
  cout << endl;
  
  // Correct MGHOT iteration count. (process starts with MGLOQD solve)
  itersMGHOT = itersMGHOT - 1;
   
  // Write iteration countrs to console 
  cout << "MGHOT iterations: " << itersMGHOT << endl;
  cout << "MGLOQD iterations: " << itersMGLOQD << endl;
  cout << "ELOT iterations: " << itersELOT << endl;
 
  // Output iteration counts 
  if (outputVars)
  {
    mesh->output->write(outputDir,"MGHOT_iters",itersMGHOT);
    mesh->output->write(outputDir,"MGLOQD_iters",itersMGLOQD);
    mesh->output->write(outputDir,"ELOT_iters",itersELOT);
    mesh->output->write(outputDir,"iterates",iters);
    mesh->output->write(outputDir,"flux_residuals",fluxResiduals);
    mesh->output->write(outputDir,"temp_residuals",tempResiduals);
    mesh->output->write(outputDir,"MGHOT_Time",mghotDuration);
    mesh->output->write(outputDir,"MGLOQD_Time",mgloqdDuration);
    mesh->output->write(outputDir,"ELOT_Time",elotDuration);
  }

  return true;

};
//==============================================================================

//==============================================================================
/// Perform a steady state solve
///
void MultilevelCoupling::solveSteadyStateResidualBalance(bool outputVars)
{

  Eigen::VectorXd xCurrentIter, xLastMGHOTIter, xLastMGLOQDIter,xLastELOTIter,\
    residualVector;
  vector<double> lastResidualELOT, lastResidualMGLOQD, residualMGHOT = {1,1,1},\
    residualMGLOQD = {1,1,1}, residualELOT = {1,1,1};
  bool eddingtonConverged;
  bool convergedMGHOT=false, convergedMGLOQD=false, convergedELOT=false;
  int itersMGHOT = 0, itersMGLOQD = 0, itersELOT = 0;
  vector<int> iters;
  vector<double> tempResMGHOT,tempResMGLOQD,tempResELOT,tempResiduals;
  vector<double> fluxResMGHOT,fluxResMGLOQD,fluxResELOT,fluxResiduals;
  vector<double> kHist;
  double power,kdiff;
 
  // Timing variables 
  double duration,totalDuration = 0.0,elotDuration = 0,\
    mgloqdDuration = 0, mghotDuration = 0;
  clock_t startTime;

  Eigen::MatrixXd volume,omega,oldFlux,newFlux;

  // Get volume and omega in each cell
  volume.setZero(mesh->nZ,mesh->nR);
  omega.setZero(mesh->nZ,mesh->nR);
  for (int iZ = 0; iZ < volume.rows(); iZ++)
  {
    for (int iR = 0; iR < volume.cols(); iR++)
    {
      volume(iZ,iR) = mesh->getGeoParams(iR,iZ)[0];
      omega(iZ,iR) = mats->omega(iZ,iR);
    }
  }
        
  // Store initial guess for k
  kHist.push_back(mats->oneGroupXS->keff);
  
  // Write mesh info 
  mesh->writeVars();
  
  while (not convergedMGHOT){ 

    ////////////////////
    // MGHOT SOLUTION // 
    ////////////////////

    // Only solve the MGHOT after we've got an estimate for the ELOT solution
    if (itersMGHOT != 0 and not p1Approx)
    {
      // Solve MGHOT problem
      cout << "MGHOT solve...";
      //startTime = clock(); 
      auto begin = chrono::high_resolution_clock::now();
      solveSteadyStateMGHOT();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      //duration = (clock() - startTime)/(double)CLOCKS_PER_SEC;
      totalDuration = totalDuration + duration; 
      mghotDuration = mghotDuration + duration; 
      cout << " done. ("<< duration << " seconds)" << endl;
      iters.push_back(3);

      // Calculate Eddington factors for MGQD problem
      cout << "Calculating MGLOQD Eddington factors...";
      eddingtonConverged = MGTToMGQD->calcEddingtonFactors();
      cout << " done."<< endl;

      // Calculate BCs for MGQD problem 
      cout << "Calculating MGLOQD boundary conditions...";
      MGTToMGQD->calcBCs();
      cout << " done."<< endl;
    }

    // Store last iterate of ELOT solution used in MGHOT level
    xLastMGHOTIter = mpqd->x;

    while (not convergedMGLOQD)
    {

      /////////////////////
      // MGLOQD SOLUTION //
      /////////////////////

      // Solve MGLOQD problem
      cout << "    ";
      cout << "MGLOQD solve..." << endl;
      //startTime = clock(); 
      auto begin = chrono::high_resolution_clock::now();
      solveSteadyStateMGLOQD();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      //duration = (clock() - startTime)/(double)CLOCKS_PER_SEC;
      totalDuration = totalDuration + duration; 
      mgloqdDuration = mgloqdDuration + duration; 
      cout << "    ";
      cout << "MGLOQD solve done. ("<< duration << " seconds)" << endl;
      iters.push_back(2);

      // Get group fluxes to use in group collapse
      mgqd->getFluxes();

      // Store last iterate of ELOT solution used in MGLOQD level
      xLastMGLOQDIter = mpqd->x;

      while (not convergedELOT)
      {

        ///////////////////
        // ELOT SOLUTION //
        ///////////////////

        // Calculate collapsed nuclear data
        cout << "        ";
        cout << "Collapsing MGLOQD data for ELOT solve...";
        MGQDToMPQD->collapseNuclearData();
        cout << " done."<< endl;

        // Store last iterate of ELOT solution used in ELOT level
        cout << "        ";
        cout << "Storing last solution...";
        xLastELOTIter = mpqd->x;
        cout << " done."<< endl;

        // Solve ELOT problem
        cout << "        ";
        cout << "ELOT solve..." << endl;
        //startTime = clock(); 
        auto begin = chrono::high_resolution_clock::now();
        solveSteadyStateELOT(mpqd->x);
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        duration = elapsed.count()*1e-9;
        //duration = (clock() - startTime)/(double)CLOCKS_PER_SEC;
        totalDuration = totalDuration + duration; 
        elotDuration = elotDuration + duration; 
        cout << "        ";
        cout << "ELOT solve done. ("<< duration << " seconds)" << endl;
        iters.push_back(1);

        // Store newest iterate 
        xCurrentIter = mpqd->x;

        lastResidualELOT = residualELOT; 

        // Calculate and print ELOT residual  
        residualELOT = MGQDToMPQD->calcResidual(xLastELOTIter,xCurrentIter);
        residualELOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastELOTIter);
        cout << "        ";
        cout << "ELOT Residual: " << residualELOT[0]; 
        cout << ", " << residualELOT[1] << endl;

        // Update iterate counters, store residuals
        itersELOT++;
        fluxResELOT.push_back(residualELOT[0]);
        tempResELOT.push_back(residualELOT[1]);

        // Store old flux
        oldFlux = mpqd->ggqd->sFlux;

        // Update variables and get new flux
        mpqd->updateSteadyStateVarsAfterConvergence(); 
        newFlux = mpqd->ggqd->sFlux;

        // Store previous eigenvalue
        mats->oneGroupXS->kold = mats->oneGroupXS->keff;

        // Calculate new eigenvalue
        mats->oneGroupXS->keff = calcK(oldFlux,newFlux,volume,\
            mats->oneGroupXS->kold);

        // Store k
        kHist.push_back(mats->oneGroupXS->keff);
        
        // Calculate power
        power = omega.cwiseProduct(mats->oneGroupXS->sigF)\
                .cwiseProduct(mpqd->ggqd->sFlux).cwiseProduct(volume).sum();

        // Scale flux to rated power
        mpqd->ggqd->sFlux = (ratedPower/power)*mpqd->ggqd->sFlux;

        // Print eigenvalue 
        cout << "        ";
        cout << "k: " << setprecision(10) <<mats->oneGroupXS->keff;
        cout << setprecision(6) <<scientific << endl;
        cout << endl;

        // Calculate difference in past and current eigenvalue
        kdiff = abs(mats->oneGroupXS->kold - mats->oneGroupXS->keff);

        // Update temperature to evaluate nuclear data at
        mats->updateTemperature(mpqd->heat->returnCurrentTemp());

        // Check if residuals are too big or if the residuals have increased
        // from the last MGLOQD residual 
        if (residualELOT[0]/lastResidualELOT[0] > resetThreshold and\
            residualELOT[1]/lastResidualELOT[1] > resetThreshold) 
        {
          // Jump back to MGLOQD level
          break;
        }


        // Check converge criteria 
        if (eps(residualMGLOQD[0], relaxTolELOT) > residualELOT[0] and\
            eps(residualMGLOQD[1], relaxTolELOT) > residualELOT[1] and\
            relaxedEpsK(residualMGLOQD[0], relaxTolELOT) > kdiff)
        {
          convergedELOT = true;
        }

        // Check keff converge criteria 
        if (abs(mats->oneGroupXS->keff - mats->oneGroupXS->kold) < 1E-10) 
        {
          //break;
        }
      } // ELOT

      // Reset convergence indicator
      convergedELOT = false; 

      // Store one group fluxes and DNPs for MGHOT and MGLOQD sources 
      mpqd->mgdnp->getCumulativeDNPDecaySource();

      lastResidualMGLOQD = residualMGLOQD;

      // Calculate and print MGLOQD residual 
      residualMGLOQD = MGQDToMPQD->calcResidual(xLastMGLOQDIter,xCurrentIter);
      residualMGLOQD = MGQDToMPQD->calcResidual(xCurrentIter,xLastMGLOQDIter);
      cout << endl;
      cout << "    ";
      cout << "MGLOQD Residual: " << residualMGLOQD[0];
      cout << ", " << residualMGLOQD[1] << endl;
      cout << endl;

      // Update iterate counters, store residuals
      itersMGLOQD++;

      fluxResMGLOQD.push_back(residualMGLOQD[0]);
      fluxResMGLOQD.insert(fluxResMGLOQD.end(),\
          fluxResELOT.begin(),fluxResELOT.end());
      fluxResELOT.clear();

      tempResMGLOQD.push_back(residualMGLOQD[1]);
      tempResMGLOQD.insert(tempResMGLOQD.end(),\
          tempResELOT.begin(),tempResELOT.end());
      tempResELOT.clear();

      // Check if residuals are too big or if the residuals have increased
      // from the last MGLOQD residual 
      if (residualMGLOQD[0]/lastResidualMGLOQD[0] > resetThreshold or\
          residualMGLOQD[1]/lastResidualMGLOQD[1] > resetThreshold)
      {
        // Jump back to MGHOT level
        break;
      }

      // Check converge criteria 
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
    residualMGHOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastMGHOTIter);
    cout << "MGHOT Residual: " << residualMGHOT[0];
    cout << ", " << residualMGHOT[1] << endl;
    cout << endl;

    // Update iterate counters, store residuals
    itersMGHOT++;

    if (itersMGHOT != 1) fluxResMGHOT.push_back(residualMGHOT[0]);
    fluxResMGHOT.insert(fluxResMGHOT.end(),\
        fluxResMGLOQD.begin(),fluxResMGLOQD.end());
    fluxResMGLOQD.clear();

    if (itersMGHOT != 1) tempResMGHOT.push_back(residualMGHOT[1]);
    tempResMGHOT.insert(tempResMGHOT.end(),\
        tempResMGLOQD.begin(),tempResMGLOQD.end());
    tempResMGLOQD.clear();

    // Update persistent iteration record
    fluxResiduals.insert(fluxResiduals.end(),fluxResMGHOT.begin(),\
        fluxResMGHOT.end());
    fluxResMGHOT.clear();

    tempResiduals.insert(tempResiduals.end(),tempResMGHOT.begin(),\
        tempResMGHOT.end());
    tempResMGHOT.clear();

    // Check converge criteria 
    if (eps(mpqd->epsMPQD) > residualMGHOT[0] and\
        eps(mpqd->epsMPQD) > residualMGHOT[1] and\
        relaxedEpsK(mpqd->epsMPQD) > kdiff) 
    {
      convergedMGHOT = true;
    }

  } //MGHOT


  // Write vars
  mpqd->updateSteadyStateVarsAfterConvergence(); 
  mgqd->updateSteadyStateVarsAfterConvergence(); 
  mpqd->writeVars(); 
  mgqd->writeVars(); 
  mats->oneGroupXS->writeVars();
  
  // Correct MGHOT iteration count. (process starts with MGLOQD solve)
  itersMGHOT = itersMGHOT - 1;

  cout << "MGHOT iterations: " << itersMGHOT << endl;
  cout << "MGLOQD iterations: " << itersMGLOQD << endl;
  cout << "ELOT iterations: " << itersELOT << endl;
 
  // Output iteration counts 
  if (outputVars)
  {
    mesh->output->write(outputDir,"Solve_Time",totalDuration);
    mesh->output->write(outputDir,"MGHOT_Time",mghotDuration);
    mesh->output->write(outputDir,"MGLOQD_Time",mgloqdDuration);
    mesh->output->write(outputDir,"ELOT_Time",elotDuration);
    mesh->output->write(outputDir,"MGHOT_iters",itersMGHOT);
    mesh->output->write(outputDir,"MGLOQD_iters",itersMGLOQD);
    mesh->output->write(outputDir,"ELOT_iters",itersELOT);
    mesh->output->write(outputDir,"iterates",iters);
    mesh->output->write(outputDir,"flux_residuals",fluxResiduals);
    mesh->output->write(outputDir,"temp_residuals",tempResiduals);
    mesh->output->write(outputDir,"k_history",kHist);
  }

};
//==============================================================================

//==============================================================================
void MultilevelCoupling::solveSteadyStateTransientResidualBalance(bool outputVars)
{
  mesh->state=0;
  solveSteadyStateResidualBalance(outputVars);
  mesh->advanceOneTimeStep();
  solveTransient();
}
//==============================================================================

//==============================================================================
/// Calculate a new eigenvalue using fission source weighting
///
double MultilevelCoupling::calcK(Eigen::MatrixXd oldFlux,\
    Eigen::MatrixXd newFlux, Eigen::MatrixXd volume, double kold)
{

  double knew;
  Eigen::MatrixXd oldFS, newFS, fsCoeff = mats->oneGroupXS->qdFluxCoeff;

  // calculate old and new fission sources weighted by new fission source
  oldFS = fsCoeff.cwiseProduct(oldFlux)\
          .cwiseProduct(fsCoeff).cwiseProduct(newFlux).cwiseProduct(volume);
  newFS = fsCoeff.cwiseProduct(newFlux)\
          .cwiseProduct(fsCoeff).cwiseProduct(newFlux).cwiseProduct(volume);
  knew = newFS.sum()/((1.0/kold)*oldFS.sum()); 

  return knew;

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
/// Dynamically determines a convergence threshold 
///
double MultilevelCoupling::relaxedEpsK(double residual, double relaxationTolerance)
{

  double eps;

  if (residual*relaxationTolerance > mpqd->epsMPQD)
    eps = residual*relaxationTolerance;
  else
    eps = epsK; 

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
/// Perform a steady state solve at the MGHOT level 
///
void MultilevelCoupling::solveSteadyStateMGHOT()
{

  // Calculate transport sources
  mgt->calcSources("steady_state");

  // Calculate transport alphas
  mgt->calcAlphas("noPrint","steady_state");

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
  if (iterativeMGLOQD)
    mgqd->solveLinearSystemIterative();
  else
    mgqd->solveLinearSystem();

  // Build neutron current system
  mgqd->buildBackCalcSystem();

  // Solve neutron current system
  mgqd->backCalculateCurrent();

};
//==============================================================================

//==============================================================================
/// Perform a steady state solve at the MGLOQD level 
///
void MultilevelCoupling::solveSteadyStateMGLOQD()
{

  // Build flux system
  mgqd->buildSteadyStateLinearSystem();

  // Solve flux system
  if (iterativeMGLOQD)
    mgqd->solveLinearSystemIterative();
  else
    mgqd->solveLinearSystem();

  // Build neutron current system
  mgqd->buildSteadyStateBackCalcSystem();

  // Solve neutron current system
  mgqd->backCalculateCurrent();

};
//==============================================================================


//==============================================================================
/// Perform a solve at the ELOT level 
///
void MultilevelCoupling::solveELOT(Eigen::VectorXd xGuess)
{

  // Build ELOT system
  mpqd->buildLinearSystem();

  // Solve ELOT system
  if (iterativeELOT)
    mpqd->solveLinearSystemIterative(xGuess);
  else
    mpqd->solveLinearSystem();
  
  //cout << "ELOT x" << endl;
  //cout << mpqd->x << endl;

};
//==============================================================================

//==============================================================================
/// Perform a steady state solve at the ELOT level 
///
void MultilevelCoupling::solveSteadyStateELOT(Eigen::VectorXd xGuess)
{

  // Build ELOT system
  mpqd->buildSteadyStateLinearSystem();

  // Solve ELOT system
  if (iterativeELOT)
    mpqd->solveLinearSystemIterative(xGuess);
  else
    mpqd->solveLinearSystem();

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
void MultilevelCoupling::solveTransient()
{

  // Timing variables
  double duration,totalDuration = 0.0;
  clock_t startTime;

  // Write mesh info
  mesh->writeVars();

  // Initialize solve 
  cout << "Computing initial solve..." << endl;
  cout << endl;
  //MGQDToMPQD->solveOneStep();
  //initialSolve();
  cout << "Initial solve completed." << endl;
  cout << endl;

  auto outerBegin = chrono::high_resolution_clock::now();

  for (int iTime = 0; iTime < mesh->dts.size(); iTime++)
  {
    cout << "Solve for t = "<< mesh->ts[iTime+1] << endl;
    cout << endl;

    //startTime = clock(); 
    auto begin = chrono::high_resolution_clock::now();
    if(solveOneStepResidualBalance(mesh->outputOnStep[iTime]))
    {
      // Report solution time
      //duration = (clock() - startTime)/(double)CLOCKS_PER_SEC;
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      cout << "Solution computed in " << duration << " seconds." << endl;      

      // Output and update variables
      mgqd->updateVarsAfterConvergence(); 
      mpqd->updateVarsAfterConvergence(); 
      if (mesh->outputOnStep[iTime])
      {
        mgqd->writeVars();
        mpqd->writeVars(); 
        mats->oneGroupXS->writeVars();
        mesh->output->write(outputDir,"Solve_Time",duration);
      }
      mesh->advanceOneTimeStep();
    }  
    else 
    {
      // Report solution time
      //duration = (clock() - startTime)/(double)CLOCKS_PER_SEC;
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      cout << "Solution aborted after " << duration << " seconds." << endl;      
      mesh->output->write(outputDir,"Solve_Time",duration);

      cout << "Solve not converged." << endl;
      cout << "Transient aborted." << endl;
      break;
    }
  } 
  auto outerEnd = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(outerEnd - outerBegin);
  duration = elapsed.count()*1e-9;
  cout << "Outer solve time: " << duration << " seconds." << endl;      

  // Report total solve time
  cout << "Total solve time: " << totalDuration << " seconds." << endl;      
  mesh->output->write(outputDir,"Solve_Time",duration,true);

};
//==============================================================================

/* PETSC FUNCTIONS */

// STEADY STATE

//==============================================================================
/// Perform a steady state solve
///
void MultilevelCoupling::solveSteadyStateResidualBalance_p(bool outputVars)
{

  Eigen::VectorXd xCurrentIter, xLastMGHOTIter, xLastMGLOQDIter,xLastELOTIter,\
    residualVector;
  vector<double> lastResidualELOT, lastResidualMGLOQD, residualMGHOT = {1,1,1},\
    residualMGLOQD = {1,1,1}, residualELOT = {1,1,1};
  bool eddingtonConverged;
  bool convergedMGHOT=false, convergedMGLOQD=false, convergedELOT=false;
  int itersMGHOT = 0, itersMGLOQD = 0, itersELOT = 0;
  vector<int> iters;
  vector<double> tempResMGHOT,tempResMGLOQD,tempResELOT,tempResiduals;
  vector<double> fluxResMGHOT,fluxResMGLOQD,fluxResELOT,fluxResiduals;
  vector<double> kHist;
  double power,kdiff;
 
  // Timing variables 
  double duration,totalDuration = 0.0,elotDuration = 0,\
    mgloqdDuration = 0, mghotDuration = 0;
  clock_t startTime;

  Eigen::MatrixXd volume,omega,oldFlux,newFlux;

  // Get volume and omega in each cell
  volume.setZero(mesh->nZ,mesh->nR);
  omega.setZero(mesh->nZ,mesh->nR);
  for (int iZ = 0; iZ < volume.rows(); iZ++)
  {
    for (int iR = 0; iR < volume.cols(); iR++)
    {
      volume(iZ,iR) = mesh->getGeoParams(iR,iZ)[0];
      omega(iZ,iR) = mats->omega(iZ,iR);
    }
  }
        
  // Store initial guess for k
  kHist.push_back(mats->oneGroupXS->keff);
  
  // Write mesh info 
  mesh->writeVars();
  
  while (not convergedMGHOT){ 

    ////////////////////
    // MGHOT SOLUTION // 
    ////////////////////

    // Only solve the MGHOT after we've got an estimate for the ELOT solution
    if (itersMGHOT != 0 and not p1Approx)
    {
      // Solve MGHOT problem
      cout << "MGHOT solve...";
      //startTime = clock(); 
      auto begin = chrono::high_resolution_clock::now();
      solveSteadyStateMGHOT();
      //duration = (clock() - startTime)/(double)CLOCKS_PER_SEC;
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      mghotDuration = mghotDuration + duration; 
      cout << " done. ("<< duration << " seconds)" << endl;
      iters.push_back(3);

      // Calculate Eddington factors for MGQD problem
      cout << "Calculating MGLOQD Eddington factors...";
      eddingtonConverged = MGTToMGQD->calcEddingtonFactors();
      cout << " done."<< endl;

      // Calculate BCs for MGQD problem 
      cout << "Calculating MGLOQD boundary conditions...";
      MGTToMGQD->calcBCs();
      cout << " done."<< endl;
    }

    // Store last iterate of ELOT solution used in MGHOT level
    //xLastMGHOTIter = mpqd->x_p;
    petscVecToEigenVec(&(mpqd->x_p),&xLastMGHOTIter);
    

    while (not convergedMGLOQD)
    {

      /////////////////////
      // MGLOQD SOLUTION //
      /////////////////////

      // Solve MGLOQD problem
      cout << "    ";
      cout << "MGLOQD solve..." << endl;
      //startTime = clock(); 
      auto begin = chrono::high_resolution_clock::now();
      solveSteadyStateMGLOQD_p();
      //duration = (clock() - startTime)/(double)CLOCKS_PER_SEC;
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      mgloqdDuration = mgloqdDuration + duration; 
      cout << "    ";
      cout << "MGLOQD solve done. ("<< duration << " seconds)" << endl;
      iters.push_back(2);

      // Get group fluxes to use in group collapse
      mgqd->getFluxes();

      // Store last iterate of ELOT solution used in MGLOQD level
      //xLastMGLOQDIter = mpqd->x_p;
      petscVecToEigenVec(&(mpqd->x_p),&xLastMGLOQDIter);

      while (not convergedELOT)
      {

        ///////////////////
        // ELOT SOLUTION //
        ///////////////////

        // Calculate collapsed nuclear data
        cout << "        ";
        cout << "Collapsing MGLOQD data for ELOT solve...";
        MGQDToMPQD->collapseNuclearData();
        cout << " done."<< endl;

        // Store last iterate of ELOT solution used in ELOT level
        cout << "        ";
        cout << "Storing last solution...";
        //xLastELOTIter = mpqd->x_p;
        petscVecToEigenVec(&(mpqd->x_p),&xLastELOTIter);
        cout << " done."<< endl;

        // Solve ELOT problem
        cout << "        ";
        cout << "ELOT solve..." << endl;
        //startTime = clock(); 
        auto begin = chrono::high_resolution_clock::now();
        solveSteadyStateELOT_p();
        //duration = (clock() - startTime)/(double)CLOCKS_PER_SEC;
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        duration = elapsed.count()*1e-9;
        totalDuration = totalDuration + duration; 
        elotDuration = elotDuration + duration; 
        cout << "        ";
        cout << "ELOT solve done. ("<< duration << " seconds)" << endl;
        iters.push_back(1);

        // Store newest iterate 
        //xCurrentIter = mpqd->x_p;
        petscVecToEigenVec(&(mpqd->x_p),&xCurrentIter);

        lastResidualELOT = residualELOT; 

        // Calculate and print ELOT residual  
        
        residualELOT = MGQDToMPQD->calcResidual(xLastELOTIter,xCurrentIter);
        residualELOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastELOTIter);
        cout << "        ";
        cout << "ELOT Residual: " << residualELOT[0]; 
        cout << ", " << residualELOT[1] << endl;

        // Update iterate counters, store residuals
        itersELOT++;
        fluxResELOT.push_back(residualELOT[0]);
        tempResELOT.push_back(residualELOT[1]);

        // Store old flux
        oldFlux = mpqd->ggqd->sFlux;

        // Update variables and get new flux
        mpqd->updateSteadyStateVarsAfterConvergence_p(); 
        newFlux = mpqd->ggqd->sFlux;

        // Store previous eigenvalue
        mats->oneGroupXS->kold = mats->oneGroupXS->keff;

        // Calculate new eigenvalue
        mats->oneGroupXS->keff = calcK(oldFlux,newFlux,volume,\
            mats->oneGroupXS->kold);

        // Store k
        kHist.push_back(mats->oneGroupXS->keff);
        
        // Calculate power
        power = omega.cwiseProduct(mats->oneGroupXS->sigF)\
                .cwiseProduct(mpqd->ggqd->sFlux).cwiseProduct(volume).sum();

        // Scale flux to rated power
        mpqd->ggqd->sFlux = (ratedPower/power)*mpqd->ggqd->sFlux;

        // Print eigenvalue 
        cout << "        ";
        cout << "k: " << setprecision(10) <<mats->oneGroupXS->keff;
        cout << setprecision(6) <<scientific << endl;
        cout << endl;

        // Calculate difference in past and current eigenvalue
        kdiff = abs(mats->oneGroupXS->kold - mats->oneGroupXS->keff);

        // Update temperature to evaluate nuclear data at
        mats->updateTemperature(mpqd->heat->returnCurrentTemp());

        // Check if residuals are too big or if the residuals have increased
        // from the last MGLOQD residual 
        if (residualELOT[0]/lastResidualELOT[0] > resetThreshold and\
            residualELOT[1]/lastResidualELOT[1] > resetThreshold) 
        {
          // Jump back to MGLOQD level
          break;
          
        }


        // Check converge criteria 
        if (eps(residualMGLOQD[0], relaxTolELOT) > residualELOT[0] and\
            eps(residualMGLOQD[1], relaxTolELOT) > residualELOT[1] and\
            relaxedEpsK(residualMGLOQD[0], relaxTolELOT) > kdiff)
        {
          convergedELOT = true;
        }

        // Check keff converge criteria 
        if (abs(mats->oneGroupXS->keff - mats->oneGroupXS->kold) < 1E-10) 
        {
          //break;
        }
      } // ELOT

      // Reset convergence indicator
      convergedELOT = false; 

      // Store one group fluxes and DNPs for MGHOT and MGLOQD sources 
      mpqd->mgdnp->getCumulativeDNPDecaySource();

      lastResidualMGLOQD = residualMGLOQD;

      // Calculate and print MGLOQD residual 
      residualMGLOQD = MGQDToMPQD->calcResidual(xLastMGLOQDIter,xCurrentIter);
      residualMGLOQD = MGQDToMPQD->calcResidual(xCurrentIter,xLastMGLOQDIter);
      cout << endl;
      cout << "    ";
      cout << "MGLOQD Residual: " << residualMGLOQD[0];
      cout << ", " << residualMGLOQD[1] << endl;
      cout << endl;

      // Update iterate counters, store residuals
      itersMGLOQD++;

      fluxResMGLOQD.push_back(residualMGLOQD[0]);
      fluxResMGLOQD.insert(fluxResMGLOQD.end(),\
          fluxResELOT.begin(),fluxResELOT.end());
      fluxResELOT.clear();

      tempResMGLOQD.push_back(residualMGLOQD[1]);
      tempResMGLOQD.insert(tempResMGLOQD.end(),\
          tempResELOT.begin(),tempResELOT.end());
      tempResELOT.clear();

      // Check if residuals are too big or if the residuals have increased
      // from the last MGLOQD residual 
      if (residualMGLOQD[0]/lastResidualMGLOQD[0] > resetThreshold or\
          residualMGLOQD[1]/lastResidualMGLOQD[1] > resetThreshold)
      {
        // Jump back to MGHOT level
        break;
      }

      // Check converge criteria 
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
    residualMGHOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastMGHOTIter);
    cout << "MGHOT Residual: " << residualMGHOT[0];
    cout << ", " << residualMGHOT[1] << endl;
    cout << endl;

    // Update iterate counters, store residuals
    itersMGHOT++;

    if (itersMGHOT != 1) fluxResMGHOT.push_back(residualMGHOT[0]);
    fluxResMGHOT.insert(fluxResMGHOT.end(),\
        fluxResMGLOQD.begin(),fluxResMGLOQD.end());
    fluxResMGLOQD.clear();

    if (itersMGHOT != 1) tempResMGHOT.push_back(residualMGHOT[1]);
    tempResMGHOT.insert(tempResMGHOT.end(),\
        tempResMGLOQD.begin(),tempResMGLOQD.end());
    tempResMGLOQD.clear();

    // Update persistent iteration record
    fluxResiduals.insert(fluxResiduals.end(),fluxResMGHOT.begin(),\
        fluxResMGHOT.end());
    fluxResMGHOT.clear();

    tempResiduals.insert(tempResiduals.end(),tempResMGHOT.begin(),\
        tempResMGHOT.end());
    tempResMGHOT.clear();

    // Check converge criteria 
    if (eps(mpqd->epsMPQD) > residualMGHOT[0] and\
        eps(mpqd->epsMPQD) > residualMGHOT[1] and\
        relaxedEpsK(mpqd->epsMPQD) > kdiff) 
    {
      convergedMGHOT = true;
    }

  } //MGHOT


  // Write vars
  mpqd->updateSteadyStateVarsAfterConvergence_p(); 
  mgqd->updateSteadyStateVarsAfterConvergence(); 
  mpqd->writeVars(); 
  mgqd->writeVars(); 
  mats->oneGroupXS->writeVars();
  
  // Correct MGHOT iteration count. (process starts with MGLOQD solve)
  itersMGHOT = itersMGHOT - 1;

  cout << "MGHOT iterations: " << itersMGHOT << endl;
  cout << "MGLOQD iterations: " << itersMGLOQD << endl;
  cout << "ELOT iterations: " << itersELOT << endl;
 
  // Output iteration counts 
  if (outputVars)
  {
    mesh->output->write(outputDir,"Solve_Time",totalDuration);
    mesh->output->write(outputDir,"MGHOT_Time",mghotDuration);
    mesh->output->write(outputDir,"MGLOQD_Time",mgloqdDuration);
    mesh->output->write(outputDir,"ELOT_Time",elotDuration);
    mesh->output->write(outputDir,"MGHOT_iters",itersMGHOT);
    mesh->output->write(outputDir,"MGLOQD_iters",itersMGLOQD);
    mesh->output->write(outputDir,"ELOT_iters",itersELOT);
    mesh->output->write(outputDir,"iterates",iters);
    mesh->output->write(outputDir,"flux_residuals",fluxResiduals);
    mesh->output->write(outputDir,"temp_residuals",tempResiduals);
    mesh->output->write(outputDir,"k_history",kHist);
  }

};
//==============================================================================

//==============================================================================
/// Perform a steady state solve at the MGLOQD level 
///
void MultilevelCoupling::solveSteadyStateMGLOQD_p()
{

  // Build flux system
  mgqd->buildSteadyStateLinearSystem_p();

  // Solve flux system
  mgqd->solveLinearSystem_p();

  // Build neutron current system
  mgqd->buildSteadyStateBackCalcSystem_p();

  // Solve neutron current system
  mgqd->backCalculateCurrent_p();

};
//==============================================================================

//==============================================================================
/// Perform a steady state solve at the ELOT level 
///
void MultilevelCoupling::solveSteadyStateELOT_p()
{

  // Build ELOT system
  mpqd->buildSteadyStateLinearSystem_p();

  // Solve ELOT system
  mpqd->solve_p();
};
//==============================================================================

// TRANSIENT

//==============================================================================
void MultilevelCoupling::solveSteadyStateTransientResidualBalance_p(bool outputVars)
{
  mesh->state=0;
  solveSteadyStateResidualBalance_p(outputVars);
  mesh->advanceOneTimeStep();
  solveTransient_p();
}
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
void MultilevelCoupling::solveTransient_p()
{

  // Timing variables
  double duration,totalDuration = 0.0;
  clock_t startTime;

  // Write mesh info
  mesh->writeVars();
  
  auto outerBegin = chrono::high_resolution_clock::now();

  for (int iTime = 0; iTime < mesh->dts.size(); iTime++)
  {
    cout << "Solve for t = "<< mesh->ts[iTime+1] << endl;
    cout << endl;

    auto begin = chrono::high_resolution_clock::now();
    if(solveOneStepResidualBalance_p(mesh->outputOnStep[iTime]))
    {
      // Report solution time
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      cout << "Solution computed in " << duration << " seconds." << endl;      

      // Output and update variables
      mgqd->updateVarsAfterConvergence(); 
      mpqd->updateVarsAfterConvergence_p(); 
      if (mesh->outputOnStep[iTime])
      {
        mgqd->writeVars();
        mpqd->writeVars(); 
        mats->oneGroupXS->writeVars();
        mesh->output->write(outputDir,"Solve_Time",duration);
      }
      mesh->advanceOneTimeStep();
    }  
    else 
    {
      // Report solution time
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      cout << "Solution aborted after " << duration << " seconds." << endl;      
      mesh->output->write(outputDir,"Solve_Time",duration);

      cout << "Solve not converged." << endl;
      cout << "Transient aborted." << endl;
      break;
    }
  } 
  
  auto outerEnd = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(outerEnd - outerBegin);
  duration = elapsed.count()*1e-9;
  cout << "Outer solve time: " << duration << " seconds." << endl;      

  // Report total solve time
  cout << "Total solve time: " << totalDuration << " seconds." << endl;      
  mesh->output->write(outputDir,"Solve_Time",duration,true);

};
//==============================================================================

//==============================================================================
/// Collapse nuclear data with flux and current weighting 
///
bool MultilevelCoupling::solveOneStepResidualBalance_p(bool outputVars)
{

  Eigen::VectorXd xCurrentIter, xLastMGHOTIter, xLastMGLOQDIter,xLastELOTIter,\
    residualVector;
  vector<double> lastResidualELOT, lastResidualMGLOQD, residualMGHOT = {1,1,1},\
    residualMGLOQD = {1,1,1}, residualELOT = {1,1,1};
  bool eddingtonConverged;
  bool convergedMGHOT=false, convergedMGLOQD=false, convergedELOT=false;
  int itersMGHOT = 0, itersMGLOQD = 0, itersELOT = 0;
  vector<int> iters;
  vector<double> tempResMGHOT,tempResMGLOQD,tempResELOT,tempResiduals;
  vector<double> fluxResMGHOT,fluxResMGLOQD,fluxResELOT,fluxResiduals;
 
  // Timing variables 
  double duration,totalDuration = 0.0,elotDuration = 0,\
    mgloqdDuration = 0, mghotDuration = 0;
  clock_t startTime;

  while (not convergedMGHOT){ 

    ////////////////////
    // MGHOT SOLUTION // 
    ////////////////////

    // Only solve the MGHOT after we've got an estimate for the ELOT solution
    if (itersMGHOT != 0 and not p1Approx)
    {
      // Solve MGHOT problem
      cout << "MGHOT solve...";
      auto begin = chrono::high_resolution_clock::now();
      solveMGHOT();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      mghotDuration = mghotDuration + duration;
      cout << " done. ("<< duration << " seconds)" << endl;
      iters.push_back(3);

      // Calculate Eddington factors for MGQD problem
      cout << "Calculating MGLOQD Eddington factors...";
      eddingtonConverged = MGTToMGQD->calcEddingtonFactors();
      cout << " done."<< endl;

      // Calculate BCs for MGQD problem 
      cout << "Calculating MGLOQD boundary conditions...";
      MGTToMGQD->calcBCs();
      cout << " done."<< endl;
    }

    // Store last iterate of ELOT solution used in MGHOT level
    petscVecToEigenVec(&(mpqd->x_p),&xLastMGHOTIter);
    
    while (not convergedMGLOQD){

      /////////////////////
      // MGLOQD SOLUTION //
      /////////////////////

      // Solve MGLOQD problem
      cout << "    ";
      cout << "MGLOQD solve..." << endl;
      auto begin = chrono::high_resolution_clock::now();
      solveMGLOQD_p();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      mgloqdDuration = mgloqdDuration + duration;
      cout << "    ";
      cout << "MGLOQD solve done. ("<< duration << " seconds)" << endl;
      iters.push_back(2);

      // Get group fluxes to use in group collapse
      mgqd->getFluxes();
      
      // Store last iterate of ELOT solution used in MGLOQD level
      petscVecToEigenVec(&(mpqd->x_p),&xLastMGLOQDIter);
    
      while (not convergedELOT) {
        
        ///////////////////
        // ELOT SOLUTION //
        ///////////////////

        // Calculate collapsed nuclear data
        cout << "        ";
        cout << "Collapsing MGLOQD data for ELOT solve...";
        MGQDToMPQD->collapseNuclearData();
        cout << " done."<< endl;
        
        // Store last iterate of ELOT solution used in ELOT level
        cout << "        ";
        cout << "Storing last solution...";
        petscVecToEigenVec(&(mpqd->x_p),&xLastELOTIter);
        cout << " done."<< endl;
      
        // Solve ELOT problem
        cout << "        ";
        cout << "ELOT solve..." << endl;
        auto begin = chrono::high_resolution_clock::now();
        solveELOT_p();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        duration = elapsed.count()*1e-9;
        elotDuration = elotDuration + duration;
        cout << "        ";
        cout << "ELOT solve done. ("<< duration << " seconds)" << endl;
        iters.push_back(1);

        // Store newest iterate 
        petscVecToEigenVec(&(mpqd->x_p),&xCurrentIter);

        lastResidualELOT = residualELOT; 

        // Calculate and print ELOT residual  
        residualELOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastELOTIter);
        cout << "        ";
        cout << "ELOT Residual: " << residualELOT[0]; 
        cout << ", " << residualELOT[1] << endl;

        // Update iterate counters, store residuals
        itersELOT++;
        fluxResELOT.push_back(residualELOT[0]);
        tempResELOT.push_back(residualELOT[1]);

        // Calculate collapsed nuclear data at new temperature
        mats->updateTemperature(mpqd->heat->returnCurrentTemp());

        // Check if residuals are too big or if the residuals have increased
        // from the last MGLOQD residual 
        if (residualELOT[0]/lastResidualELOT[0] > resetThreshold and\
            residualELOT[1]/lastResidualELOT[1] > resetThreshold) 
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
      residualMGLOQD = MGQDToMPQD->calcResidual(xCurrentIter,xLastMGLOQDIter);
      cout << endl;
      cout << "    ";
      cout << "MGLOQD Residual: " << residualMGLOQD[0];
      cout << ", " << residualMGLOQD[1] << endl;
      cout << endl;
        
      // Update iterate counters, store residuals
      itersMGLOQD++;

      fluxResMGLOQD.push_back(residualMGLOQD[0]);
      fluxResMGLOQD.insert(fluxResMGLOQD.end(),\
          fluxResELOT.begin(),fluxResELOT.end());
      fluxResELOT.clear();

      tempResMGLOQD.push_back(residualMGLOQD[1]);
      tempResMGLOQD.insert(tempResMGLOQD.end(),\
          tempResELOT.begin(),tempResELOT.end());
      tempResELOT.clear();

      // Check if residuals are too big or if the residuals have increased
      // from the last MGLOQD residual 
      if (residualMGLOQD[0]/lastResidualMGLOQD[0] > resetThreshold or\
          residualMGLOQD[1]/lastResidualMGLOQD[1] > resetThreshold) 
      {
        // Jump back to MGHOT level
        break;
      }

      // Check converge criteria 
      if (eps(residualMGHOT[0],relaxTolMGLOQD) > residualMGLOQD[0] and\
          eps(residualMGHOT[1],relaxTolMGLOQD) > residualMGLOQD[1])
      { 
        convergedMGLOQD = true;
      }

    } // MGLOQD
    
    // Reset convergence indicator
    convergedMGLOQD = false;
     
    // Calculate and print MGHOT residual 
    residualMGHOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastMGHOTIter);
    cout << "MGHOT Residual: " << residualMGHOT[0];
    cout << ", " << residualMGHOT[1] << endl;
    cout << endl;
      
    // Update iterate counters, store residuals
    itersMGHOT++;
    
    if (itersMGHOT != 1) fluxResMGHOT.push_back(residualMGHOT[0]);
    fluxResMGHOT.insert(fluxResMGHOT.end(),\
        fluxResMGLOQD.begin(),fluxResMGLOQD.end());
    fluxResMGLOQD.clear();
    
    if (itersMGHOT != 1) tempResMGHOT.push_back(residualMGHOT[1]);
    tempResMGHOT.insert(tempResMGHOT.end(),\
        tempResMGLOQD.begin(),tempResMGLOQD.end());
    tempResMGLOQD.clear();

    // Update persistent iteration record
    fluxResiduals.insert(fluxResiduals.end(),fluxResMGHOT.begin(),\
        fluxResMGHOT.end());
    fluxResMGHOT.clear();
    
    tempResiduals.insert(tempResiduals.end(),tempResMGHOT.begin(),\
        tempResMGHOT.end());
    tempResMGHOT.clear();
      
    // Check converge criteria 
    if (eps(mpqd->epsMPQD) > residualMGHOT[0] and\
        eps(mpqd->epsMPQD) > residualMGHOT[1]) 
    {
      convergedMGHOT = true;
    }

  } //MGHOT
    
  cout << endl;

  // Correct MGHOT iteration count. (process starts with MGLOQD solve)
  itersMGHOT = itersMGHOT - 1;
   
  // Write iteration countrs to console 
  cout << "MGHOT iterations: " << itersMGHOT << endl;
  cout << "MGLOQD iterations: " << itersMGLOQD << endl;
  cout << "ELOT iterations: " << itersELOT << endl;
 
  // Output iteration counts 
  if (outputVars)
  {
    mesh->output->write(outputDir,"MGHOT_iters",itersMGHOT);
    mesh->output->write(outputDir,"MGLOQD_iters",itersMGLOQD);
    mesh->output->write(outputDir,"ELOT_iters",itersELOT);
    mesh->output->write(outputDir,"iterates",iters);
    mesh->output->write(outputDir,"flux_residuals",fluxResiduals);
    mesh->output->write(outputDir,"temp_residuals",tempResiduals);
    mesh->output->write(outputDir,"MGHOT_Time",mghotDuration);
    mesh->output->write(outputDir,"MGLOQD_Time",mgloqdDuration);
    mesh->output->write(outputDir,"ELOT_Time",elotDuration);
  }

  return true;

};
//==============================================================================

//==============================================================================
/// Perform a solve at the MGLOQD level 
///
void MultilevelCoupling::solveMGLOQD_p()
{
  
  PetscErrorCode ierr;

  // Build flux system
  mgqd->buildLinearSystem_p();

  // Solve flux system
  mgqd->solveLinearSystem_p();
  
  // Build neutron current system
  mgqd->buildBackCalcSystem_p();
  
  //cout << "MGLOQD b_p" << endl;
  //VecView(mgqd->QDSolve->b_p,PETSC_VIEWER_STDOUT_WORLD);

  // Solve neutron current system
  mgqd->backCalculateCurrent_p();

};
//==============================================================================

//==============================================================================
/// Perform a solve at the ELOT level 
///
void MultilevelCoupling::solveELOT_p()
{

  // Build ELOT system
  mpqd->buildLinearSystem_p();

  mpqd->solve_p();
  
  //cout << "ELOT x_p" << endl;
  //VecView(mpqd->x_p,PETSC_VIEWER_STDOUT_WORLD);

};
//==============================================================================

// PSUEDO TRANSIENT

//==============================================================================
void MultilevelCoupling::solveSteadyStatePsuedoTransient_p(bool outputVars)
{
  mesh->state=0;
  solveSteadyStateResidualBalance_p(outputVars);
  mesh->advanceOneTimeStep();
  solvePsuedoTransient_p();
}
//==============================================================================

//==============================================================================
/// Solve multiple step psuedo transient
///
void MultilevelCoupling::solvePsuedoTransient_p()
{

  // Timing variables
  double duration,totalDuration = 0.0;
  clock_t startTime;

  // Write mesh info
  mesh->writeVars();

  // Initialize solve 

  auto outerBegin = chrono::high_resolution_clock::now();

  for (int iTime = 0; iTime < mesh->dts.size(); iTime++)
  {
    cout << "Solve for t = "<< mesh->ts[iTime+1] << endl;
    cout << endl;

    //startTime = clock(); 
    auto begin = chrono::high_resolution_clock::now();
    if(solvePsuedoTransientResidualBalance_p(mesh->outputOnStep[iTime]))
    {
      // Report solution time
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      cout << "Solution computed in " << duration << " seconds." << endl;      

      // Output and update variables
      //mgqd->updateVarsAfterConvergence(); 
      //mpqd->updateVarsAfterConvergence(); 
      if (mesh->outputOnStep[iTime])
      {
        mgqd->writeVars();
        mpqd->writeVars(); 
        mats->oneGroupXS->writeVars();
        mesh->output->write(outputDir,"Solve_Time",duration);
      }
      mesh->advanceOneTimeStep();
    }  
    else 
    {
      // Report solution time
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      cout << "Solution aborted after " << duration << " seconds." << endl;      
      mesh->output->write(outputDir,"Solve_Time",duration);

      cout << "Solve not converged." << endl;
      cout << "Transient aborted." << endl;
      break;
    }
  } 
  auto outerEnd = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(outerEnd - outerBegin);
  duration = elapsed.count()*1e-9;
  cout << "Outer solve time: " << duration << " seconds." << endl;      

  // Report total solve time
  cout << "Total solve time: " << totalDuration << " seconds." << endl;      
  mesh->output->write(outputDir,"Solve_Time",duration,true);

};
//==============================================================================

//==============================================================================
/// Perform psuedo transient solve for a single step
///
bool MultilevelCoupling::solvePsuedoTransientResidualBalance_p(bool outputVars)
{

  Eigen::VectorXd xCurrentIter, xLastMGHOTIter, xLastMGLOQDIter,xLastELOTIter,\
    residualVector;
  vector<double> lastResidualELOT, lastResidualMGLOQD, residualMGHOT = {1,1,1},\
    residualMGLOQD = {1,1,1}, residualELOT = {1,1,1};
  bool eddingtonConverged;
  bool convergedMGHOT=false, convergedMGLOQD=false, convergedELOT=false;
  int itersMGHOT = 0, itersMGLOQD = 0, itersELOT = 0;
  vector<int> iters;
  vector<double> tempResMGHOT,tempResMGLOQD,tempResELOT,tempResiduals;
  vector<double> fluxResMGHOT,fluxResMGLOQD,fluxResELOT,fluxResiduals;
  vector<double> kHist;
  double power,kdiff;
 
  // Timing variables 
  double duration,totalDuration = 0.0,elotDuration = 0,\
    mgloqdDuration = 0, mghotDuration = 0;
  clock_t startTime;

  Eigen::MatrixXd volume,omega,oldFlux,newFlux;

  // Get volume and omega in each cell
  volume.setZero(mesh->nZ,mesh->nR);
  omega.setZero(mesh->nZ,mesh->nR);
  for (int iZ = 0; iZ < volume.rows(); iZ++)
  {
    for (int iR = 0; iR < volume.cols(); iR++)
    {
      volume(iZ,iR) = mesh->getGeoParams(iR,iZ)[0];
      omega(iZ,iR) = mats->omega(iZ,iR);
    }
  }
        
  // Store initial guess for k
  mats->oneGroupXS->keff=1.0;
  kHist.push_back(mats->oneGroupXS->keff);
  
  // Write mesh info 
  mesh->writeVars();
  
  while (not convergedMGHOT){ 

    ////////////////////
    // MGHOT SOLUTION // 
    ////////////////////

    // Only solve the MGHOT after we've got an estimate for the ELOT solution
    if (itersMGHOT != 0 and not p1Approx)
    {
      // Solve MGHOT problem
      cout << "MGHOT solve...";
      auto begin = chrono::high_resolution_clock::now();
      solveSteadyStateMGHOT();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      mghotDuration = mghotDuration + duration; 
      cout << " done. ("<< duration << " seconds)" << endl;
      iters.push_back(3);

      // Calculate Eddington factors for MGQD problem
      cout << "Calculating MGLOQD Eddington factors...";
      eddingtonConverged = MGTToMGQD->calcEddingtonFactors();
      cout << " done."<< endl;

      // Calculate BCs for MGQD problem 
      cout << "Calculating MGLOQD boundary conditions...";
      MGTToMGQD->calcBCs();
      cout << " done."<< endl;
    }

    // Store last iterate of ELOT solution used in MGHOT level
    petscVecToEigenVec(&(mpqd->x_p),&xLastMGHOTIter);

    while (not convergedMGLOQD)
    {

      /////////////////////
      // MGLOQD SOLUTION //
      /////////////////////

      // Solve MGLOQD problem
      cout << "    ";
      cout << "MGLOQD solve..." << endl;
      auto begin = chrono::high_resolution_clock::now();
      solveSteadyStateMGLOQD_p();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      mgloqdDuration = mgloqdDuration + duration; 
      cout << "    ";
      cout << "MGLOQD solve done. ("<< duration << " seconds)" << endl;
      iters.push_back(2);

      // Get group fluxes to use in group collapse
      mgqd->getFluxes();

      // Store last iterate of ELOT solution used in MGLOQD level
      petscVecToEigenVec(&(mpqd->x_p),&xLastMGLOQDIter);

      while (not convergedELOT)
      {

        ///////////////////
        // ELOT SOLUTION //
        ///////////////////

        // Calculate collapsed nuclear data
        cout << "        ";
        cout << "Collapsing MGLOQD data for ELOT solve...";
        MGQDToMPQD->collapseNuclearData();
        cout << " done."<< endl;

        // Store last iterate of ELOT solution used in ELOT level
        cout << "        ";
        cout << "Storing last solution...";
        //xLastELOTIter = mpqd->x_p;
        petscVecToEigenVec(&(mpqd->x_p),&xLastELOTIter);
        cout << " done."<< endl;

        // Solve ELOT problem
        cout << "        ";
        cout << "ELOT solve..." << endl;
        auto begin = chrono::high_resolution_clock::now();
        solvePsuedoTransientELOT_p();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        duration = elapsed.count()*1e-9;
        totalDuration = totalDuration + duration; 
        elotDuration = elotDuration + duration; 
        cout << "        ";
        cout << "ELOT solve done. ("<< duration << " seconds)" << endl;
        iters.push_back(1);

        // Store newest iterate 
        petscVecToEigenVec(&(mpqd->x_p),&xCurrentIter);

        lastResidualELOT = residualELOT; 

        // Calculate and print ELOT residual  
        
        residualELOT = MGQDToMPQD->calcResidual(xLastELOTIter,xCurrentIter);
        residualELOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastELOTIter);
        cout << "        ";
        cout << "ELOT Residual: " << residualELOT[0]; 
        cout << ", " << residualELOT[1] << endl;

        // Update iterate counters, store residuals
        itersELOT++;
        fluxResELOT.push_back(residualELOT[0]);
        tempResELOT.push_back(residualELOT[1]);

        // Store old flux
        oldFlux = mpqd->ggqd->sFlux;

        // Update variables and get new flux
        mpqd->updateSteadyStateVarsAfterConvergence_p(); 
        newFlux = mpqd->ggqd->sFlux;

        // Store previous eigenvalue
        mats->oneGroupXS->kold = mats->oneGroupXS->keff;

        // Calculate new eigenvalue
        mats->oneGroupXS->keff = calcK(oldFlux,newFlux,volume,\
            mats->oneGroupXS->kold);

        // Store k
        kHist.push_back(mats->oneGroupXS->keff);
        
        // Calculate power
        power = omega.cwiseProduct(mats->oneGroupXS->sigF)\
                .cwiseProduct(mpqd->ggqd->sFlux).cwiseProduct(volume).sum();

        // Scale flux to rated power
        mpqd->ggqd->sFlux = (ratedPower/power)*mpqd->ggqd->sFlux;

        // Print eigenvalue 
        cout << "        ";
        cout << "k: " << setprecision(10) <<mats->oneGroupXS->keff;
        cout << setprecision(6) <<scientific << endl;
        cout << endl;

        // Calculate difference in past and current eigenvalue
        kdiff = abs(mats->oneGroupXS->kold - mats->oneGroupXS->keff);

        // Update temperature to evaluate nuclear data at
        mats->updateTemperature(mpqd->heat->returnCurrentTemp());

        // Check if residuals are too big or if the residuals have increased
        // from the last MGLOQD residual 
        if (residualELOT[0]/lastResidualELOT[0] > resetThreshold and\
            residualELOT[1]/lastResidualELOT[1] > resetThreshold) 
        {
          // Jump back to MGLOQD level
          break;
          
        }


        // Check converge criteria 
        if (eps(residualMGLOQD[0], relaxTolELOT) > residualELOT[0] and\
            eps(residualMGLOQD[1], relaxTolELOT) > residualELOT[1] and\
            relaxedEpsK(residualMGLOQD[0], relaxTolELOT) > kdiff)
        {
          convergedELOT = true;
        }

        // Check keff converge criteria 
        if (abs(mats->oneGroupXS->keff - mats->oneGroupXS->kold) < 1E-10) 
        {
          //break;
        }
      } // ELOT

      // Reset convergence indicator
      convergedELOT = false; 

      // Store one group fluxes and DNPs for MGHOT and MGLOQD sources 
      mpqd->mgdnp->getCumulativeDNPDecaySource();

      lastResidualMGLOQD = residualMGLOQD;

      // Calculate and print MGLOQD residual 
      residualMGLOQD = MGQDToMPQD->calcResidual(xLastMGLOQDIter,xCurrentIter);
      residualMGLOQD = MGQDToMPQD->calcResidual(xCurrentIter,xLastMGLOQDIter);
      cout << endl;
      cout << "    ";
      cout << "MGLOQD Residual: " << residualMGLOQD[0];
      cout << ", " << residualMGLOQD[1] << endl;
      cout << endl;

      // Update iterate counters, store residuals
      itersMGLOQD++;

      fluxResMGLOQD.push_back(residualMGLOQD[0]);
      fluxResMGLOQD.insert(fluxResMGLOQD.end(),\
          fluxResELOT.begin(),fluxResELOT.end());
      fluxResELOT.clear();

      tempResMGLOQD.push_back(residualMGLOQD[1]);
      tempResMGLOQD.insert(tempResMGLOQD.end(),\
          tempResELOT.begin(),tempResELOT.end());
      tempResELOT.clear();

      // Check if residuals are too big or if the residuals have increased
      // from the last MGLOQD residual 
      if (residualMGLOQD[0]/lastResidualMGLOQD[0] > resetThreshold or\
          residualMGLOQD[1]/lastResidualMGLOQD[1] > resetThreshold)
      {
        // Jump back to MGHOT level
        break;
      }

      // Check converge criteria 
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
    residualMGHOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastMGHOTIter);
    cout << "MGHOT Residual: " << residualMGHOT[0];
    cout << ", " << residualMGHOT[1] << endl;
    cout << endl;

    // Update iterate counters, store residuals
    itersMGHOT++;

    if (itersMGHOT != 1) fluxResMGHOT.push_back(residualMGHOT[0]);
    fluxResMGHOT.insert(fluxResMGHOT.end(),\
        fluxResMGLOQD.begin(),fluxResMGLOQD.end());
    fluxResMGLOQD.clear();

    if (itersMGHOT != 1) tempResMGHOT.push_back(residualMGHOT[1]);
    tempResMGHOT.insert(tempResMGHOT.end(),\
        tempResMGLOQD.begin(),tempResMGLOQD.end());
    tempResMGLOQD.clear();

    // Update persistent iteration record
    fluxResiduals.insert(fluxResiduals.end(),fluxResMGHOT.begin(),\
        fluxResMGHOT.end());
    fluxResMGHOT.clear();

    tempResiduals.insert(tempResiduals.end(),tempResMGHOT.begin(),\
        tempResMGHOT.end());
    tempResMGHOT.clear();

    // Check converge criteria 
    if (eps(mpqd->epsMPQD) > residualMGHOT[0] and\
        eps(mpqd->epsMPQD) > residualMGHOT[1] and\
        relaxedEpsK(mpqd->epsMPQD) > kdiff) 
    {
      convergedMGHOT = true;
    }

  } //MGHOT


  // Write vars
  mpqd->updateSteadyStateVarsAfterConvergence_p(); 
  mgqd->updateSteadyStateVarsAfterConvergence(); 
  mpqd->writeVars(); 
  mgqd->writeVars(); 
  mats->oneGroupXS->writeVars();
  
  // Correct MGHOT iteration count. (process starts with MGLOQD solve)
  itersMGHOT = itersMGHOT - 1;

  cout << "MGHOT iterations: " << itersMGHOT << endl;
  cout << "MGLOQD iterations: " << itersMGLOQD << endl;
  cout << "ELOT iterations: " << itersELOT << endl;
 
  // Output iteration counts 
  if (outputVars)
  {
    mesh->output->write(outputDir,"Solve_Time",totalDuration);
    mesh->output->write(outputDir,"MGHOT_Time",mghotDuration);
    mesh->output->write(outputDir,"MGLOQD_Time",mgloqdDuration);
    mesh->output->write(outputDir,"ELOT_Time",elotDuration);
    mesh->output->write(outputDir,"MGHOT_iters",itersMGHOT);
    mesh->output->write(outputDir,"MGLOQD_iters",itersMGLOQD);
    mesh->output->write(outputDir,"ELOT_iters",itersELOT);
    mesh->output->write(outputDir,"iterates",iters);
    mesh->output->write(outputDir,"flux_residuals",fluxResiduals);
    mesh->output->write(outputDir,"temp_residuals",tempResiduals);
    mesh->output->write(outputDir,"k_history",kHist);
  }

  return true;

};
//==============================================================================

//==============================================================================
/// Perform a solve at the ELOT level 
///
void MultilevelCoupling::solvePsuedoTransientELOT_p()
{

  // Build ELOT system
  mpqd->buildPsuedoTransientLinearSystem_p();

  mpqd->solve_p();
  
  //cout << "ELOT x_p" << endl;
  //VecView(mpqd->x_p,PETSC_VIEWER_STDOUT_WORLD);

};
//==============================================================================

//==============================================================================
/// Read in optional parameters that might be specified in the input 
///
void MultilevelCoupling::checkOptionalParameters()
{
  string boundaryType;

  // Check for relaxTolELOT specification..
  if ((*input)["parameters"]["relaxTolELOT"])
    relaxTolELOT=(*input)["parameters"]["relaxTolELOT"].as<double>();

  // Check for relaxTolMGLOQD specification.
  if ((*input)["parameters"]["relaxTolELOT"])
    relaxTolMGLOQD=(*input)["parameters"]["relaxTolMGLOQD"].as<double>();

  // Check for resetThreshold specification.
  if ((*input)["parameters"]["resetThreshold"])
    resetThreshold=(*input)["parameters"]["resetThreshold"].as<double>();

  // Check for rated power 
  if ((*input)["parameters"]["ratedPower"])
    ratedPower=(*input)["parameters"]["ratedPower"].as<double>();

  // Check for k convergence criteria
  if ((*input)["parameters"]["epsK"])
    epsK=(*input)["parameters"]["epsK"].as<double>();

  // Check if iterative solver should be used for ELOT 
  if ((*input)["parameters"]["iterativeELOT"])
    iterativeELOT=(*input)["parameters"]["iterativeELOT"].as<bool>();

  // Check if iterative solver should be used for MGLOQD 
  if ((*input)["parameters"]["iterativeMGLOQD"])
    iterativeMGLOQD=(*input)["parameters"]["iterativeMGLOQD"].as<bool>();

  // Check if the P1 approximation should be used
  if ((*input)["parameters"]["mgqd-bcs"])
  {
    boundaryType=(*input)["parameters"]["mgqd-bcs"].as<string>();

    if (boundaryType == "diffusion" or boundaryType == "DIFFUSION" \
        or boundaryType == "Diffusion")
      p1Approx = true;
  } 


};
//==============================================================================
