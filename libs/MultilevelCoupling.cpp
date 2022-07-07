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
double MultilevelCoupling::epsTemp(double residual, double relaxationTolerance)
{

  double eps;

  if (residual*relaxationTolerance > mpqd->epsMPQDTemp)
    eps = residual*relaxationTolerance;
  else
    eps = mpqd->epsMPQDTemp; 

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

/* PETSC FUNCTIONS */

// STEADY STATE

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
  
  auto outerBegin = chrono::high_resolution_clock::now();
  
  while (not convergedMGHOT){ 

    ////////////////////
    // MGHOT SOLUTION // 
    ////////////////////

    // Only solve the MGHOT after we've got an estimate for the ELOT solution
    if (itersMGHOT != 0 and not p1Approx)
    {
      // Solve MGHOT problem
      PetscPrintf(PETSC_COMM_WORLD,"MGHOT solve...");
      auto begin = chrono::high_resolution_clock::now();
      solveSteadyStateMGHOT();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      mghotDuration = mghotDuration + duration; 
      PetscPrintf(PETSC_COMM_WORLD," done. (%6.4e seconds)\n",duration);
      iters.push_back(3);

      // Calculate Eddington factors for MGQD problem
      PetscPrintf(PETSC_COMM_WORLD,"Calculating MGLOQD Eddington factors...");
      eddingtonConverged = MGTToMGQD->calcEddingtonFactors();
      PetscPrintf(PETSC_COMM_WORLD," done.\n");

      // Calculate BCs for MGQD problem 
      PetscPrintf(PETSC_COMM_WORLD,"Calculating MGLOQD boundary conditions...");
      MGTToMGQD->calcBCs();
      PetscPrintf(PETSC_COMM_WORLD," done.\n\n");
    }

    // Store last iterate of ELOT solution used in MGHOT level
    petscVecToEigenVec(&(mpqd->x_p),&xLastMGHOTIter);
    

    while (not convergedMGLOQD)
    {

      /////////////////////
      // MGLOQD SOLUTION //
      /////////////////////

      // Solve MGLOQD problem
      PetscPrintf(PETSC_COMM_WORLD,"    ");
      PetscPrintf(PETSC_COMM_WORLD,"MGLOQD solve...\n");
      auto begin = chrono::high_resolution_clock::now();
      solveSteadyStateMGLOQD();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      mgloqdDuration = mgloqdDuration + duration; 
      PetscPrintf(PETSC_COMM_WORLD,"    ");
      PetscPrintf(PETSC_COMM_WORLD,"MGLOQD solve done. (%6.4e seconds)\n\n",duration);
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
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"Collapsing MGLOQD data for ELOT solve...");
        MGQDToMPQD->collapseNuclearData();
        PetscPrintf(PETSC_COMM_WORLD," done.\n");

        // Store last iterate of ELOT solution used in ELOT level
        petscVecToEigenVec(&(mpqd->x_p),&xLastELOTIter);

        // Solve ELOT problem
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"ELOT solve...\n");
        auto begin = chrono::high_resolution_clock::now();
        solveSteadyStateELOT();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        duration = elapsed.count()*1e-9;
        totalDuration = totalDuration + duration; 
        elotDuration = elotDuration + duration; 
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"ELOT solve done. (%6.4e seconds)\n",duration);
        iters.push_back(1);

        // Store newest iterate 
        petscVecToEigenVec(&(mpqd->x_p),&xCurrentIter);

        lastResidualELOT = residualELOT; 

        // Calculate and print ELOT residual  
        residualELOT = MGQDToMPQD->calcResidual(xLastELOTIter,xCurrentIter);
        residualELOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastELOTIter);
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"ELOT Residual: %6.4e, %6.4e \n",residualELOT[0],\
            residualELOT[1]);

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

        if (ratedPower > 0) 
          // Scale flux to rated power
          mpqd->ggqd->sFlux = (ratedPower/power)*mpqd->ggqd->sFlux;
        else
          // Normalize flux
          mpqd->ggqd->sFlux = (fluxNormalization/mpqd->ggqd->sFlux.sum())*mpqd->ggqd->sFlux;

        // Print eigenvalue 
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"k: %11.10e \n\n",mats->oneGroupXS->keff);

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
      //PetscPrintf(PETSC_COMM_WORLD,"\n");
      PetscPrintf(PETSC_COMM_WORLD,"    ");
      PetscPrintf(PETSC_COMM_WORLD,"MGLOQD Residual: %6.4e, %6.4e \n\n",residualMGLOQD[0],\
          residualMGLOQD[1]);

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
    PetscPrintf(PETSC_COMM_WORLD,"MGHOT Residual: %6.4e, %6.4e \n\n",residualMGHOT[0],\
        residualMGHOT[1]);

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

  auto outerEnd = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(outerEnd - outerBegin);
  duration = elapsed.count()*1e-9;
  //cout << "Outer solve time: " << duration << " seconds." << endl;      
  PetscPrintf(PETSC_COMM_WORLD,"Outer solve time: %6.4e\n", duration);

  // Write vars
  mpqd->updateSteadyStateVarsAfterConvergence(); 
  mgqd->updateSteadyStateVarsAfterConvergence(); 
  
  if (mesh->verbose_keff_only)
    mesh->output->write(mpqd->outputDir,"keff",mats->oneGroupXS->keff);
  else
  {
    mpqd->writeVars(); 
    mgqd->writeVars(); 
    mats->oneGroupXS->writeVars();
  }

  // Correct MGHOT iteration count. (process starts with MGLOQD solve)
  itersMGHOT = itersMGHOT - 1;

  PetscPrintf(PETSC_COMM_WORLD,"MGHOT iterations: %i\n", itersMGHOT);
  PetscPrintf(PETSC_COMM_WORLD,"MGLOQD iterations: %i\n", itersMGLOQD);
  PetscPrintf(PETSC_COMM_WORLD,"ELOT iterations: %i\n", itersELOT);

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
    mesh->output->write(mpqd->outputDir,"Power",power);
  }

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
  mgqd->solveLinearSystem();

  // Build neutron current system
  mgqd->buildSteadyStateBackCalcSystem();

  // Solve neutron current system
  mgqd->backCalculateCurrent();

};
//==============================================================================

//==============================================================================
/// Perform a steady state solve at the ELOT level 
///
void MultilevelCoupling::solveSteadyStateELOT()
{

  // Build ELOT system
  mpqd->buildSteadyStateLinearSystem();

  // Solve ELOT system
  mpqd->solve();
};
//==============================================================================

// TRANSIENT

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
/// Collapse nuclear data with flux and current weighting 
///
void MultilevelCoupling::solveTransient()
{

  // Timing variables
  double duration,totalDuration = 0.0;
  clock_t startTime;

  // Write mesh info
  mesh->writeVars();

  auto outerBegin = chrono::high_resolution_clock::now();

  for (int iTime = 0; iTime < mesh->dts.size(); iTime++)
  {
    //cout << "Solve for t = "<< mesh->ts[iTime+1] << endl;
    //cout << endl;
    PetscPrintf(PETSC_COMM_WORLD,"Solve for t = %6.4e\n\n",mesh->ts[iTime+1]);

    auto begin = chrono::high_resolution_clock::now();
    if(solveOneStepResidualBalance(mesh->outputOnStep[iTime]))
    {
      // Report solution time
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      //cout << "Solution computed in " << duration << " seconds." << endl;      
      PetscPrintf(PETSC_COMM_WORLD,"Solution computed in %6.4e seconds.\n", duration); 

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
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      //cout << "Solution aborted after " << duration << " seconds." << endl;      
      PetscPrintf(PETSC_COMM_WORLD,"Solution aborted after %6.4e seconds.\n", duration); 
      mesh->output->write(outputDir,"Solve_Time",duration);

      //cout << "Solve not converged." << endl;
      //cout << "Transient aborted." << endl;
      PetscPrintf(PETSC_COMM_WORLD,"Solve not converged\n"); 
      PetscPrintf(PETSC_COMM_WORLD,"Transient aborted.\n"); 
      break;
    }
  } 

  auto outerEnd = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(outerEnd - outerBegin);
  duration = elapsed.count()*1e-9;
  PetscPrintf(PETSC_COMM_WORLD,"Outer solve time: %6.4e\n", duration);

  // Report total solve time
  PetscPrintf(PETSC_COMM_WORLD,"Total solve time: %6.4e\n", totalDuration);
  mesh->output->write(outputDir,"Solve_Time",duration,true);

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

  while (not convergedMGHOT){ 

    ////////////////////
    // MGHOT SOLUTION // 
    ////////////////////

    // Only solve the MGHOT after we've got an estimate for the ELOT solution
    if (itersMGHOT != 0 and not p1Approx)
    {
      // Solve MGHOT problem
      PetscPrintf(PETSC_COMM_WORLD,"MGHOT solve...");
      auto begin = chrono::high_resolution_clock::now();
      solveMGHOT();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      mghotDuration = mghotDuration + duration;
      PetscPrintf(PETSC_COMM_WORLD," done. (%6.4e seconds)\n",duration);
      iters.push_back(3);

      // Calculate Eddington factors for MGQD problem
      PetscPrintf(PETSC_COMM_WORLD,"Calculating MGLOQD Eddington factors...");
      eddingtonConverged = MGTToMGQD->calcEddingtonFactors();
      PetscPrintf(PETSC_COMM_WORLD," done.\n");

      // Calculate BCs for MGQD problem 
      PetscPrintf(PETSC_COMM_WORLD,"Calculating MGLOQD boundary conditions...");
      MGTToMGQD->calcBCs();
      PetscPrintf(PETSC_COMM_WORLD," done.\n\n");
    }

    // Store last iterate of ELOT solution used in MGHOT level
    petscVecToEigenVec(&(mpqd->x_p),&xLastMGHOTIter);

    while (not convergedMGLOQD){

      /////////////////////
      // MGLOQD SOLUTION //
      /////////////////////

      // Solve MGLOQD problem
      PetscPrintf(PETSC_COMM_WORLD,"    ");
      PetscPrintf(PETSC_COMM_WORLD,"MGLOQD solve...\n");
      auto begin = chrono::high_resolution_clock::now();
      solveMGLOQD();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      mgloqdDuration = mgloqdDuration + duration;
      PetscPrintf(PETSC_COMM_WORLD,"    ");
      PetscPrintf(PETSC_COMM_WORLD,"MGLOQD solve done. (%6.4e seconds)\n\n",duration);
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
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"Collapsing MGLOQD data for ELOT solve...");
        MGQDToMPQD->collapseNuclearData();
        PetscPrintf(PETSC_COMM_WORLD," done.\n");

        // Store last iterate of ELOT solution used in ELOT level
        petscVecToEigenVec(&(mpqd->x_p),&xLastELOTIter);

        // Solve ELOT problem
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"ELOT solve...\n");
        auto begin = chrono::high_resolution_clock::now();
        solveELOT();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        duration = elapsed.count()*1e-9;
        elotDuration = elotDuration + duration;
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"ELOT solve done. (%6.4e seconds)\n",duration);
        iters.push_back(1);

        // Store newest iterate 
        petscVecToEigenVec(&(mpqd->x_p),&xCurrentIter);

        lastResidualELOT = residualELOT; 

        // Calculate and print ELOT residual  
        residualELOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastELOTIter);
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"ELOT Residual: %6.4e, %6.4e \n\n",residualELOT[0],\
            residualELOT[1]);

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
      PetscPrintf(PETSC_COMM_WORLD,"    ");
      PetscPrintf(PETSC_COMM_WORLD,"MGLOQD Residual: %6.4e, %6.4e \n\n",residualMGLOQD[0],\
          residualMGLOQD[1]);

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
    PetscPrintf(PETSC_COMM_WORLD,"MGHOT Residual: %6.4e, %6.4e \n\n",residualMGHOT[0],\
        residualMGHOT[1]);

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
  PetscPrintf(PETSC_COMM_WORLD,"MGHOT iterations: %i\n", itersMGHOT);
  PetscPrintf(PETSC_COMM_WORLD,"MGLOQD iterations: %i\n", itersMGLOQD);
  PetscPrintf(PETSC_COMM_WORLD,"ELOT iterations: %i\n", itersELOT);

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
void MultilevelCoupling::solveMGLOQD()
{
  PetscErrorCode ierr;

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

  mpqd->solve();
};
//==============================================================================

// PSEUDO TRANSIENT

//==============================================================================
void MultilevelCoupling::solveSteadyStatePseudoTransient(bool outputVars)
{
  mesh->state=0;
  MatSetOption(mpqd->A_p, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  solveSteadyStateResidualBalance(outputVars);
  mesh->advanceOneTimeStep();
  MatDestroy(&(mpqd->A_p));
  initPETScMat(&(mpqd->A_p),mpqd->nUnknowns,40);
  solvePseudoTransient();
}
//==============================================================================

//==============================================================================
/// Solve multiple step pseudo transient
///
void MultilevelCoupling::solvePseudoTransient()
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
    //cout << "Solve for t = "<< mesh->ts[iTime+1] << endl;
    //cout << endl;
    PetscPrintf(PETSC_COMM_WORLD,"Solve for t = %6.4e\n\n",mesh->ts[iTime+1]);

    // Evaluate material velocites at current time step
    mats->readFlowVelocity(mesh->ts[iTime+1]);

    //startTime = clock(); 
    auto begin = chrono::high_resolution_clock::now();
    if(solvePseudoTransientResidualBalance(mesh->outputOnStep[iTime]))
    {
      // Report solution time
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      //cout << "Solution computed in " << duration << " seconds." << endl;      
      PetscPrintf(PETSC_COMM_WORLD,"Solution computed in %6.4e seconds.\n", duration); 

      // Output and update variables
      if (mesh->outputOnStep[iTime])
      {

        if (mesh->verbose_keff_only)
        {
          mesh->output->write(mpqd->outputDir,"keff",mats->oneGroupXS->keff);
        }
        else
        {
          mgqd->writeVars();
          mpqd->writeVars(); 
          mats->oneGroupXS->writeVars();
        }
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
      //cout << "Solution aborted after " << duration << " seconds." << endl;      
      PetscPrintf(PETSC_COMM_WORLD,"Solution aborted after %6.4e seconds.\n", duration); 
      mesh->output->write(outputDir,"Solve_Time",duration);

      //cout << "Solve not converged." << endl;
      //cout << "Transient aborted." << endl;
      PetscPrintf(PETSC_COMM_WORLD,"Solve not converged\n"); 
      PetscPrintf(PETSC_COMM_WORLD,"Transient aborted.\n"); 
      break;
    }
  } 
  auto outerEnd = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(outerEnd - outerBegin);
  duration = elapsed.count()*1e-9;
  //cout << "Outer solve time: " << duration << " seconds." << endl;      
  PetscPrintf(PETSC_COMM_WORLD,"Outer solve time: %6.4e\n", duration);

  // Report total solve time
  //cout << "Total solve time: " << totalDuration << " seconds." << endl;      
  PetscPrintf(PETSC_COMM_WORLD,"Total solve time: %6.4e\n", totalDuration);
  mesh->output->write(outputDir,"Solve_Time",duration,true);

};
//==============================================================================

//==============================================================================
/// Perform pseudo transient solve for a single step
///
bool MultilevelCoupling::solvePseudoTransientResidualBalance(bool outputVars)
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
      PetscPrintf(PETSC_COMM_WORLD,"MGHOT solve...");
      auto begin = chrono::high_resolution_clock::now();
      solveSteadyStateMGHOT();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      mghotDuration = mghotDuration + duration; 
      PetscPrintf(PETSC_COMM_WORLD," done. (%6.4e seconds)\n",duration);
      iters.push_back(3);

      // Calculate Eddington factors for MGQD problem
      PetscPrintf(PETSC_COMM_WORLD,"Calculating MGLOQD Eddington factors...");
      eddingtonConverged = MGTToMGQD->calcEddingtonFactors();
      PetscPrintf(PETSC_COMM_WORLD," done.\n");

      // Calculate BCs for MGQD problem 
      PetscPrintf(PETSC_COMM_WORLD,"Calculating MGLOQD boundary conditions...");
      MGTToMGQD->calcBCs();
      PetscPrintf(PETSC_COMM_WORLD," done.\n\n");
    }

    // Store last iterate of ELOT solution used in MGHOT level
    petscVecToEigenVec(&(mpqd->x_p),&xLastMGHOTIter);

    while (not convergedMGLOQD)
    {

      /////////////////////
      // MGLOQD SOLUTION //
      /////////////////////

      // Solve MGLOQD problem
      PetscPrintf(PETSC_COMM_WORLD,"    ");
      PetscPrintf(PETSC_COMM_WORLD,"MGLOQD solve...\n");
      auto begin = chrono::high_resolution_clock::now();
      solveSteadyStateMGLOQD();
      auto end = chrono::high_resolution_clock::now();
      auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
      duration = elapsed.count()*1e-9;
      totalDuration = totalDuration + duration; 
      mgloqdDuration = mgloqdDuration + duration; 
      PetscPrintf(PETSC_COMM_WORLD,"    ");
      PetscPrintf(PETSC_COMM_WORLD,"MGLOQD solve done. (%6.4e seconds)\n\n",duration);
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
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"Collapsing MGLOQD data for ELOT solve...");
        MGQDToMPQD->collapseNuclearData();
        PetscPrintf(PETSC_COMM_WORLD," done.\n");

        // Store last iterate of ELOT solution used in ELOT level
        petscVecToEigenVec(&(mpqd->x_p),&xLastELOTIter);

        // Solve ELOT problem
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"ELOT solve...\n");
        auto begin = chrono::high_resolution_clock::now();
        solvePseudoTransientELOT();
        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        duration = elapsed.count()*1e-9;
        totalDuration = totalDuration + duration; 
        elotDuration = elotDuration + duration; 
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"ELOT solve done. (%6.4e seconds)\n",duration);
        iters.push_back(1);

        // Store newest iterate 
        petscVecToEigenVec(&(mpqd->x_p),&xCurrentIter);

        lastResidualELOT = residualELOT; 

        // Calculate and print ELOT residual  

        residualELOT = MGQDToMPQD->calcResidual(xLastELOTIter,xCurrentIter);
        residualELOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastELOTIter);
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"ELOT Residual: %6.4e, %6.4e \n\n",residualELOT[0],\
            residualELOT[1]);

        // Update iterate counters, store residuals
        itersELOT++;
        fluxResELOT.push_back(residualELOT[0]);
        tempResELOT.push_back(residualELOT[1]);

        // Store old flux
        oldFlux = mpqd->ggqd->sFlux;

        // Update variables and get new flux
        mpqd->updatePseudoTransientVars(); 
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

        if (ratedPower > 0) 
          // Scale flux to rated power
          mpqd->ggqd->sFlux = (ratedPower/power)*mpqd->ggqd->sFlux;
        else
          // Normalize flux
          mpqd->ggqd->sFlux = (fluxNormalization/mpqd->ggqd->sFlux.sum())*mpqd->ggqd->sFlux;

        // Print eigenvalue 
        PetscPrintf(PETSC_COMM_WORLD,"        ");
        PetscPrintf(PETSC_COMM_WORLD,"k: %11.10e \n\n",mats->oneGroupXS->keff);

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
        else if (not mesh->verbose)
        { 
          mesh->output->deleteLines(7);
        }


        // Check converge criteria 
        if (eps(residualMGLOQD[0], relaxTolELOT) > residualELOT[0] and\
            epsTemp(residualMGLOQD[1], relaxTolELOT) > residualELOT[1] and\
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
      PetscPrintf(PETSC_COMM_WORLD,"    ");
      PetscPrintf(PETSC_COMM_WORLD,"MGLOQD Residual: %6.4e, %6.4e \n\n",residualMGLOQD[0],\
          residualMGLOQD[1]);

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
          epsTemp(residualMGHOT[1],relaxTolMGLOQD) > residualMGLOQD[1])
      { 
        convergedMGLOQD = true;
      }

    } // MGLOQD

    // Reset convergence indicator
    convergedMGLOQD = false;

    // Calculate and print MGHOT residual 
    residualMGHOT = MGQDToMPQD->calcResidual(xLastMGHOTIter,xCurrentIter);
    residualMGHOT = MGQDToMPQD->calcResidual(xCurrentIter,xLastMGHOTIter);
    PetscPrintf(PETSC_COMM_WORLD,"MGHOT Residual: %6.4e, %6.4e \n\n",residualMGHOT[0],\
        residualMGHOT[1]);

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
        epsTemp(mpqd->epsMPQD) > residualMGHOT[1] and\
        relaxedEpsK(mpqd->epsMPQD) > kdiff) 
    {
      convergedMGHOT = true;
    }

  } //MGHOT


  // Write vars
  mpqd->updateSteadyStateVarsAfterConvergence(); 
  mgqd->updateSteadyStateVarsAfterConvergence(); 

  // Correct MGHOT iteration count. (process starts with MGLOQD solve)
  itersMGHOT = itersMGHOT - 1;

  PetscPrintf(PETSC_COMM_WORLD,"MGHOT iterations: %i\n", itersMGHOT);
  PetscPrintf(PETSC_COMM_WORLD,"MGLOQD iterations: %i\n", itersMGLOQD);
  PetscPrintf(PETSC_COMM_WORLD,"ELOT iterations: %i\n", itersELOT);

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
    mesh->output->write(mpqd->outputDir,"Power",power);
  }

  return true;

};
//==============================================================================

//==============================================================================
/// Perform a solve at the ELOT level 
///
void MultilevelCoupling::solvePseudoTransientELOT()
{

  // Build ELOT system
  mpqd->buildPseudoTransientLinearSystem();

  mpqd->solve();

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
  else
    ratedPower=-1.0;
  
  // Check for flux normalization factor
  if ((*input)["parameters"]["fluxNormalization"])
    fluxNormalization=(*input)["parameters"]["fluxNormalization"].as<double>();
  else
    fluxNormalization=1.0;

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
