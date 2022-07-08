//// File: MultiPhysicsCoupledQD.cpp     
// Purpose: Contains precursor, heat, and grey group quasidiffusion objects. 
//   Also owns and is responsible for solving the coupled system each of 
//   of those object contribute to.
// Date: April 9, 2020

#include "MultiPhysicsCoupledQD.h"
#include "HeatTransfer.h"
#include "MultiGroupDNP.h"
#include "GreyGroupQD.h"

using namespace std;

//==============================================================================
/// MultiPhysicsCoupledQD class object constructor
///
/// @param [in] blankType blank for this material
MultiPhysicsCoupledQD::MultiPhysicsCoupledQD(Materials * myMats,\
    Mesh * myMesh,\
    YAML::Node * myInput) 
{
  int tempIndexOffset;

  // Assign inputs to their member variables
  mats = myMats;
  mesh = myMesh;
  input = myInput;

  // Initialize grey group qd object
  ggqd = new GreyGroupQD(mats,mesh,input,this);
  ggqd->indexOffset = 0;  

  // Initialize heat transfer object and set index offset
  heat = new HeatTransfer(mats,mesh,input,this);
  heat->indexOffset = ggqd->nUnknowns;

  // Calculate index offset for MGDNP object initialization
  tempIndexOffset = ggqd->nUnknowns + heat->nUnknowns;
  mgdnp = new MultiGroupDNP(mats,mesh,input,this,tempIndexOffset);

  // Set size of linear system
  nUnknowns = ggqd->nUnknowns + heat->nUnknowns + mgdnp->nCoreUnknowns;

  /* Initialize PETSc variables */
  // Multiphysics system variables 
  initPETScMat(&A_p,nUnknowns,40);
  initPETScVec(&x_p,nUnknowns);
  initPETScVec(&xPast_p,nUnknowns);
  initPETScVec(&b_p,nUnknowns);

  // Initialize sequential variables
  initPETScVec(&xPast_p_seq,nUnknowns);

  // Assign pointers in ggqd object
  ggqd->GGSolver->assignMPQDPointer(this);

  // Initialize xPast 
  setInitialCondition();

  // Check optional parameters
  checkOptionalParams();
};
//==============================================================================

//==============================================================================
/// Include a flux source in the linear system
///
/// @param [in] iZ axial location 
/// @param [in] iR radial location
/// @param [in] iEq equation index
/// @param [in] coeff coefficient of flux source
int MultiPhysicsCoupledQD::fluxSource(int iZ,int iR,int iEq,double coeff)
{
  int iCF = 0; // index of cell-average flux value in index vector  
  vector<int> indices = ggqd->GGSolver->getIndices(iR,iZ);
  PetscErrorCode ierr;

  ierr = MatSetValue(A_p,iEq,indices[0],coeff,ADD_VALUES);CHKERRQ(ierr); 

  return ierr;
};
//==============================================================================

//==============================================================================
/// Include a dnp source in the linear system
///
/// @param [in] iZ axial location 
/// @param [in] iR radial location
/// @param [in] iEq equation index
/// @param [in] coeff coefficient of dnp source
int MultiPhysicsCoupledQD::dnpSource(int iZ,int iR,int iEq,double coeff)
{
  int index,indexOffset;
  double groupLambda;
  PetscErrorCode ierr;
  PetscScalar value;

  for (int iGroup = 0; iGroup < mgdnp->DNPs.size(); ++iGroup)
  {
    indexOffset = mgdnp->DNPs[iGroup]->coreIndexOffset;
    index = mgdnp->DNPs[iGroup]->getIndex(iZ,iR,indexOffset);
    groupLambda = mgdnp->DNPs[iGroup]->lambda;

    value = coeff*groupLambda;
    ierr = MatSetValue(A_p,iEq,index,value,ADD_VALUES);CHKERRQ(ierr); 
  }

  return ierr;
};
//==============================================================================

//==============================================================================
/// Map values in multiphysics objects into xPast
///
void MultiPhysicsCoupledQD::setInitialCondition()
{

  // Object for broadcasting PETSc variable 
  VecScatter     ctx;

  // Set fluxes 
  ggqd->GGSolver->setFlux();

  // Set temperatures
  heat->setTemp();

  // Set DNP concentrations in core and recirculation loop
  mgdnp->setCoreDNPConc();  
  mgdnp->setRecircDNPConc();  

  /* Calculate currents consistent with fluxes in xPast */
  x_p = xPast_p; // getFlux() pulls from x

  // Broadcast xPast
  VecScatterCreateToAll(xPast_p,&ctx,&(xPast_p_seq));
  VecScatterBegin(ctx,xPast_p,xPast_p_seq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,xPast_p,xPast_p_seq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);

  ggqd->GGSolver->getFlux();
  ggqd->GGSolver->formBackCalcSystem();
  ggqd->GGSolver->backCalculateCurrent();
  ggqd->GGSolver->getCurrent();

  VecZeroEntries(x_p);
};
//==============================================================================

//==============================================================================
/// Write variables out to CSVs
///
void MultiPhysicsCoupledQD::writeVars()
{
  string name; 

  // Scalar flux
  mesh->output->write(outputDir,"Flux",ggqd->sFlux);

  // Face fluxes
  mesh->output->write(outputDir,"Flux_Radial",ggqd->sFluxR);
  mesh->output->write(outputDir,"Flux_Axial",ggqd->sFluxZ);

  // Currents
  mesh->output->write(outputDir,"Current_Radial",ggqd->currentR);
  mesh->output->write(outputDir,"Current_Axial",ggqd->currentZ);

  // Eddington factors
  mesh->output->write(outputDir,"Err",ggqd->Err);
  mesh->output->write(outputDir,"Ezz",ggqd->Ezz);
  mesh->output->write(outputDir,"Erz",ggqd->Erz);
  mesh->output->write(outputDir,"ErrAxial",ggqd->ErrAxial);
  mesh->output->write(outputDir,"EzzAxial",ggqd->EzzAxial);
  mesh->output->write(outputDir,"ErzAxial",ggqd->ErzAxial);
  mesh->output->write(outputDir,"ErrRadial",ggqd->ErrRadial);
  mesh->output->write(outputDir,"EzzRadial",ggqd->EzzRadial);
  mesh->output->write(outputDir,"ErzRadial",ggqd->ErzRadial);
  mesh->output->write(outputDir,"GL",ggqd->GL);
  mesh->output->write(outputDir,"GR",ggqd->GR);

  // Eigenvalue
  mesh->output->write(outputDir,"keff",mats->oneGroupXS->keff);

  // DNP concentrations
  for (int iDNP = 0; iDNP < mgdnp->DNPs.size(); iDNP++)
  {
    name =  "Core_DNP_Concentration_Group_"\
             + to_string(mgdnp->DNPs[iDNP]->dnpID);
    mesh->output->write(outputDir,name,mgdnp->DNPs[iDNP]->dnpConc); 
    name =  "Recirculation_DNP_Concentration_Group_"\
             + to_string(mgdnp->DNPs[iDNP]->dnpID);
    mesh->output->write(outputDir,name,mgdnp->DNPs[iDNP]->recircConc); 
  }

  // Temperature
  mesh->output->write(outputDir,"Temperature",heat->temp);
};
//==============================================================================

//==============================================================================
/// Run transient with multiple solves 
///
void MultiPhysicsCoupledQD::printVars()
{

  // Read solutions from 1D vector to 2D matrices 
  cout << "Flux:" << endl; 
  cout << ggqd->sFlux << endl; 

  cout << "Temperature:" << endl; 
  cout << heat->temp << endl; 

  cout << "Core DNP concentration:" << endl; 
  mgdnp->printCoreDNPConc();

  cout << "Recirc DNP concentration:" << endl; 
  mgdnp->printRecircDNPConc();

};
//==============================================================================

// Dual purpose

//==============================================================================
/// Solve linear system for multiphysics coupled quasidiffusion system with a 
/// direct solve
///
int MultiPhysicsCoupledQD::solve()
{

  string PetscSolver = "bicg";
  string PetscPreconditioner = "bjacobi";
  PetscErrorCode ierr;
  int its,m,n;
  double norm,relTol,absTol;

  auto begin = chrono::high_resolution_clock::now();
  /* Get matrix dimensions */
  ierr = MatGetSize(A_p, &m, &n); CHKERRQ(ierr);
  
  // Set solve tolerances
  absTol= 1e-50;
  relTol = 1.e-9/((m+1)*(n+1));
  if (relTol < 0.0)
    relTol = 1.e-14;

  /* Create solver */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A_p,A_p);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,relTol,absTol,PETSC_DEFAULT,PETSC_DEFAULT);
  CHKERRQ(ierr);

  /* Set solver type */
  ierr = KSPSetType(ksp,PetscSolver.c_str());CHKERRQ(ierr);

  /* Set preconditioner type */
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PetscPreconditioner.c_str());CHKERRQ(ierr);

  /* Solve the system */
  ierr = KSPSolve(ksp,b_p,x_p);CHKERRQ(ierr);
  auto end = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

  /* Print solve information */
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);

  /* Destroy solver */
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  /* Solve for recirculation concentrations */
  mgdnp->solveRecircLinearSystem();

  return ierr;
};
//==============================================================================

// Steady state

//==============================================================================
/// Build linear system for multiphysics coupled quasidiffusion system
///
int MultiPhysicsCoupledQD::buildSteadyStateLinearSystem()
{
  PetscErrorCode ierr;

  // Reset linear system
  MatZeroEntries(A_p);
  VecZeroEntries(b_p);

  // Build QD system
  ggqd->buildSteadyStateLinearSystem();

  // Build heat transfer system
  heat->buildSteadyStateLinearSystem();

  // Build delayed neutron precursor balance system in core
  mgdnp->buildSteadyStateCoreLinearSystem();  

  // Build delayed neutron precursor balance system in recirculation loop
  mgdnp->buildSteadyStateRecircLinearSystem();  

  /* Finalize assembly for A_p and b_p */
  ierr = MatAssemblyBegin(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  
  ierr = VecAssemblyBegin(b_p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b_p);CHKERRQ(ierr);

  return ierr;
};
//==============================================================================

//==============================================================================
/// Update 2D variables after a converged steady state solve
///
void MultiPhysicsCoupledQD::updateSteadyStateVarsAfterConvergence()
{
  // Object for broadcasting PETSc variable 
  VecScatter     ctx;

  // Read solutions from 1D vector to 2D matrices 
  ggqd->GGSolver->getFlux();

  heat->getTemp();

  mgdnp->getCoreDNPConc();

  mgdnp->getRecircDNPConc();

  // Back calculate currents
  // ToDo Add PETSc support for back calc system
  ggqd->GGSolver->formSteadyStateBackCalcSystem();
  ggqd->GGSolver->backCalculateCurrent();
  ggqd->GGSolver->getCurrent();

  // Set xPast and past neutron velocities 
  xPast_p = x_p;
  mats->oneGroupXS->neutVPast = mats->oneGroupXS->neutV;  
  mats->oneGroupXS->zNeutVPast = mats->oneGroupXS->zNeutV;  
  mats->oneGroupXS->rNeutVPast = mats->oneGroupXS->rNeutV;  

  // Broadcast xPast
  VecDestroy(&(xPast_p_seq));
  VecScatterCreateToAll(xPast_p,&ctx,&(xPast_p_seq));
  VecScatterBegin(ctx,xPast_p,xPast_p_seq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,xPast_p,xPast_p_seq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);
};
//==============================================================================

//==============================================================================
/// Update 2D vars after a converged pseudo transient solve
///
void MultiPhysicsCoupledQD::updatePseudoTransientVars()
{
  // Object for broadcasting PETSc variable 
  VecScatter     ctx;

  // Read solutions from 1D vector to 2D matrices 
  ggqd->GGSolver->getFlux();

  // Back calculate currents
  ggqd->GGSolver->formSteadyStateBackCalcSystem();
  ggqd->GGSolver->backCalculateCurrent();
  ggqd->GGSolver->getCurrent();
};
//==============================================================================


//==============================================================================
/// Run an ELOT-only steady state 
/// Note: not accessible to user
void MultiPhysicsCoupledQD::solveSteadyState()
{
  for (int iT = 0; iT < mesh->dts.size(); iT++)
  {
    buildSteadyStateLinearSystem();
    solve();   
    updateSteadyStateVarsAfterConvergence();
  }

  writeVars();
};
//==============================================================================

/* TRANSIENT */

//==============================================================================
/// Run ELOT-only transient with multiple solves 
/// Note: not accessible to user
void MultiPhysicsCoupledQD::solveTransient()
{
  for (int iT = 0; iT < mesh->dts.size(); iT++)
  {
    buildLinearSystem();
    solve();   
    updateVarsAfterConvergence();
  }

  writeVars();

};
//==============================================================================

//==============================================================================
/// Build linear system for multiphysics coupled quasidiffusion system
///
int MultiPhysicsCoupledQD::buildLinearSystem()
{
  PetscErrorCode ierr;

  // Reset linear system
  MatZeroEntries(A_p);
  VecZeroEntries(b_p);

  // Build QD system
  ggqd->buildLinearSystem();

  // Build heat transfer system
  heat->buildLinearSystem();

  // Build delayed neutron precursor balance system in core
  mgdnp->buildCoreLinearSystem();  

  // Build delayed neutron precursor balance system in recirculation loop
  mgdnp->buildRecircLinearSystem();  

  /* Finalize assembly for A_p and b_p */
  ierr = MatAssemblyBegin(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  
  ierr = VecAssemblyBegin(b_p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b_p);CHKERRQ(ierr);

  return ierr;
};
//==============================================================================

// Build pseudo-transient system

//==============================================================================
/// Build linear system for multiphysics coupled quasidiffusion system
///
int MultiPhysicsCoupledQD::buildPseudoTransientLinearSystem()
{
  PetscErrorCode ierr;

  // Reset linear system
  MatZeroEntries(A_p);
  VecZeroEntries(b_p);

  // Build QD system
  ggqd->buildSteadyStateLinearSystem();

  // Build heat transfer system
  heat->buildLinearSystem();

  // Build delayed neutron precursor balance system in core
  mgdnp->buildPseudoTransientCoreLinearSystem();  

  // Build delayed neutron precursor balance system in recirculation loop
  mgdnp->buildRecircLinearSystem();  

  /* Finalize assembly for A_p and b_p */
  ierr = MatAssemblyBegin(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_p,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  
  ierr = VecAssemblyBegin(b_p);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b_p);CHKERRQ(ierr);

  return ierr;
};
//==============================================================================

//==============================================================================
/// Update variables after convergence
///
void MultiPhysicsCoupledQD::updateVarsAfterConvergence()
{
  // Object for broadcasting PETSc variable 
  VecScatter     ctx;

  // Read solutions from 1D vector to 2D matrices 
  ggqd->GGSolver->getFlux();

  heat->getTemp();

  mgdnp->getCoreDNPConc();

  mgdnp->getRecircDNPConc();

  // Back calculate currents
  ggqd->GGSolver->formBackCalcSystem();
  ggqd->GGSolver->backCalculateCurrent();
  ggqd->GGSolver->getCurrent();

  // Set xPast and past neutron velocities 
  xPast_p = x_p;
  mats->oneGroupXS->neutVPast = mats->oneGroupXS->neutV;  
  mats->oneGroupXS->zNeutVPast = mats->oneGroupXS->zNeutV;  
  mats->oneGroupXS->rNeutVPast = mats->oneGroupXS->rNeutV;  

  // Broadcast xPast
  VecDestroy(&(xPast_p_seq));
  VecScatterCreateToAll(xPast_p,&ctx,&(xPast_p_seq));
  VecScatterBegin(ctx,xPast_p,xPast_p_seq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,xPast_p,xPast_p_seq,\
      INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);

};
//==============================================================================

//==============================================================================
/// Check for optional input parameters of relevance to this object
void MultiPhysicsCoupledQD::checkOptionalParams()
{
  string precondInput;

  if ((*input)["parameters"]["epsFlux"])
    epsFlux = (*input)["parameters"]["epsFlux"].as<double>();

  if ((*input)["parameters"]["epsTemp"])
    epsTemp = (*input)["parameters"]["epsTemp"].as<double>();

}
//==============================================================================
