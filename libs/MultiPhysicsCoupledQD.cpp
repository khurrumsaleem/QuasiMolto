// File: MultiPhysicsCoupledQD.cpp     
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
  int tempIndexOffset,nUnknowns;

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
  A.resize(nUnknowns,nUnknowns); 
  x.setZero(nUnknowns); 
  xPast.setOnes(nUnknowns); 
  b.setZero(nUnknowns);

  // Assign pointers in ggqd object
  ggqd->GGSolver->assignPointers(&A,&x,&xPast,&b);
 
   // Initialize xPast 
  initializeXPast();

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
void MultiPhysicsCoupledQD::fluxSource(int iZ,int iR,int iEq,double coeff,\
    Eigen::SparseMatrix<double,Eigen::RowMajor> * myA)
{

  int iCF = 0; // index of cell-average flux value in index vector  
  vector<int> indices = ggqd->GGSolver->getIndices(iR,iZ);

  myA->coeffRef(iEq,indices[0]) += coeff; 

};
//==============================================================================

//==============================================================================
/// Include a flux source in the linear system
///
/// @param [in] iZ axial location 
/// @param [in] iR radial location
/// @param [in] iEq equation index
/// @param [in] coeff coefficient of flux source
void MultiPhysicsCoupledQD::fluxSource(int iZ,int iR,int iEq,double coeff,\
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> * myA)
{

  int iCF = 0; // index of cell-average flux value in index vector  
  vector<int> indices = ggqd->GGSolver->getIndices(iR,iZ);

  (*myA)(iEq,indices[0]) += coeff; 

};
//==============================================================================


//==============================================================================
/// Include a dnp source in the linear system
///
/// @param [in] iZ axial location 
/// @param [in] iR radial location
/// @param [in] iEq equation index
/// @param [in] coeff coefficient of dnp source
void MultiPhysicsCoupledQD::dnpSource(int iZ,int iR,int iEq,double coeff,\
    Eigen::SparseMatrix<double,Eigen::RowMajor> * myA)
{
  int index,indexOffset;
  double groupLambda;

  for (int iGroup = 0; iGroup < mgdnp->DNPs.size(); ++iGroup)
  {
    indexOffset = mgdnp->DNPs[iGroup]->coreIndexOffset;
    index = mgdnp->DNPs[iGroup]->getIndex(iZ,iR,indexOffset);
    groupLambda = mgdnp->DNPs[iGroup]->lambda;

    myA->coeffRef(iEq,index) += coeff*groupLambda;
  }

};
//==============================================================================

//==============================================================================
/// Build linear system for multiphysics coupled quasidiffusion system
///
void MultiPhysicsCoupledQD::buildLinearSystem()
{

  // Get number of non-zero elements in sparse matrix to optimize building 
  // linear system
  int nonZeros = A.nonZeros(); 

  // Reset linear system
  A.setZero();
  A.reserve(nonZeros); 
  x.setZero();
  b.setZero();

  // Build QD system
  ggqd->buildLinearSystem();

  // Build heat transfer system
  heat->buildLinearSystem();

  // Build delayed neutron precursor balance system in core
  mgdnp->buildCoreLinearSystem();  

  // Build delayed neutron precursor balance system in recirculation loop
  mgdnp->buildRecircLinearSystem();  

};
//==============================================================================

//==============================================================================
/// Build linear system for multiphysics coupled quasidiffusion system
///
void MultiPhysicsCoupledQD::buildSteadyStateLinearSystem()
{

  // Get number of non-zero elements in sparse matrix to optimize building 
  // linear system
  int nonZeros = A.nonZeros(); 

  // Reset linear system
  A.setZero();
  A.reserve(nonZeros); 
  x.setZero();
  b.setZero();

  // Build QD system
  ggqd->buildSteadyStateLinearSystem();

  // Build heat transfer system
  heat->buildSteadyStateLinearSystem();

  // Build delayed neutron precursor balance system in core
  mgdnp->buildSteadyStateCoreLinearSystem();  

  // Build delayed neutron precursor balance system in recirculation loop
  mgdnp->buildSteadyStateRecircLinearSystem();  

};
//==============================================================================


//==============================================================================
/// Map values in multiphysics objects into xPast
///
void MultiPhysicsCoupledQD::initializeXPast()
{

  // Set fluxes 
  ggqd->GGSolver->setFlux();

  // Set temperatures
  heat->setTemp();

  // Set DNP concentrations in core and recirculation loop
  mgdnp->setCoreDNPConc();  
  mgdnp->setRecircDNPConc();  

  // Calculate currents consistent with fluxes in xPast
  x = xPast; // getFlux() pulls from x
  ggqd->GGSolver->getFlux();
  ggqd->GGSolver->formBackCalcSystem();
  ggqd->GGSolver->backCalculateCurrent();
  ggqd->GGSolver->getCurrent();
  x.setZero(); // Reset x
};
//==============================================================================

//==============================================================================
/// Solve linear system for multiphysics coupled quasidiffusion system
///
void MultiPhysicsCoupledQD::solveLinearSystem()
{
  
  solveSuperLU();
  mgdnp->solveRecircLinearSystem();

};
//==============================================================================

//==============================================================================
/// Solve linear system for multiphysics coupled quasidiffusion system with an
/// iterative solver
///
void MultiPhysicsCoupledQD::solveLinearSystemIterative(Eigen::VectorXd xGuess)
{

  double duration,totalDuration = 0.0;
  clock_t startTime;
  int n = Eigen::nbThreads();
  int solveOutcome;
  //cout << "number procs: " << n << endl;

  if (preconditioner == iluPreconditioner) 
    solveOutcome = solveIterativeILU(xGuess);
  else if (preconditioner == diagPreconditioner)
  {
    solveOutcome = solveIterativeDiag(xGuess);
    if (solveOutcome != Eigen::Success)
    {
      cout << "            " << "BiCGSTAB solve failed! Attempting iterative";
      cout << " solve with ILU preconditioner." << endl;
      solveOutcome = solveIterativeILU(xGuess);
    }
  }

  if (solveOutcome != Eigen::Success)
  {
    cout << "            " << "Iterative solve failed! ";
    cout << "Using SuperLU direct solve.";  
    solveSuperLU();
  }

  // Solve recirc system
  mgdnp->solveRecircLinearSystem();

};
//==============================================================================

//==============================================================================
/// Solve linear system for multiphysics coupled quasidiffusion system with a 
/// direct solve
///
int MultiPhysicsCoupledQD::solveSuperLU()
{

  int success;

  // Declare SuperLU solver
  Eigen::SuperLU<Eigen::SparseMatrix<double>> solverLU;
  A.makeCompressed();
  solverLU.compute(A);
  x = solverLU.solve(b);

  // Return outcome of solve
  return success = solverLU.info();
  
};
//==============================================================================

//==============================================================================
/// Solve linear system for multiphysics coupled quasidiffusion system with an
/// iterative solver and incomplete LU preconditioner
///
int MultiPhysicsCoupledQD::solveIterativeILU(Eigen::VectorXd xGuess)
{

  int success;

  // Declare solver with ILUT preconditioner
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>,\
    Eigen::IncompleteLUT<double> > solver;

  // Set preconditioner parameters
  solver.preconditioner().setDroptol(1E-4);
  //solver.preconditioner().setFillfactor(100);
  //solver.preconditioner().setFillfactor(5);

  // Set convergence parameters
  //solver.setTolerance(1E-14);
  //solver.setMaxIterations(20);

  // Solve system
  A.makeCompressed();
  solver.analyzePattern(A);
  solver.factorize(A);
  x = solver.solveWithGuess(b,xGuess);

  if (mesh->verbose) 
  {
    cout << "            ";
    cout << "info:     " << solver.info() << endl;
    cout << "            ";
    cout << "#iterations:     " << solver.iterations() << endl;
    cout << "            ";
    cout << "estimated error: " << solver.error() << endl;
    cout << "            ";
    cout << "tolerance: " << solver.tolerance() << endl;
  }

  // Return outcome of solve
  return success = solver.info();
  
};
//==============================================================================

//==============================================================================
/// Solve linear system for multiphysics coupled quasidiffusion system with an
/// iterative solver and diagonal preconditioner
///
int MultiPhysicsCoupledQD::solveIterativeDiag(Eigen::VectorXd xGuess)
{

  int success;

  // Declare solver with default diagonal precondition (cheaper to calculate) 
  // but usually requires more iterations to converge
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor> > solver;

  // Set convergence parameters
  //solver.setTolerance(1E-14);
  //solver.setMaxIterations(20);

  // Solve system
  A.makeCompressed();
  solver.analyzePattern(A);
  solver.factorize(A);
  x = solver.solveWithGuess(b,xGuess);

  if (mesh->verbose) 
  {
    cout << "            ";
    cout << "info:     " << solver.info() << endl;
    cout << "            ";
    cout << "#iterations:     " << solver.iterations() << endl;
    cout << "            ";
    cout << "estimated error: " << solver.error() << endl;
    cout << "            ";
    cout << "tolerance: " << solver.tolerance() << endl;
  }

  // Return outcome of solve
  return success = solver.info();

};
//==============================================================================

//==============================================================================
/// Run transient with multiple solves 
///
void MultiPhysicsCoupledQD::updateVarsAfterConvergence()
{

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
  xPast = x;
  mats->oneGroupXS->neutVPast = mats->oneGroupXS->neutV;  
  mats->oneGroupXS->zNeutVPast = mats->oneGroupXS->zNeutV;  
  mats->oneGroupXS->rNeutVPast = mats->oneGroupXS->rNeutV;  

};
//==============================================================================

//==============================================================================
/// Run transient with multiple solves 
///
void MultiPhysicsCoupledQD::updateSteadyStateVarsAfterConvergence()
{

  // Read solutions from 1D vector to 2D matrices 
  ggqd->GGSolver->getFlux();

  heat->getTemp();

  mgdnp->getCoreDNPConc();

  mgdnp->getRecircDNPConc();

  // Back calculate currents
  ggqd->GGSolver->formSteadyStateBackCalcSystem();
  ggqd->GGSolver->backCalculateCurrent();
  ggqd->GGSolver->getCurrent();

  // Set xPast and past neutron velocities 
  xPast = x;
  mats->oneGroupXS->neutVPast = mats->oneGroupXS->neutV;  
  mats->oneGroupXS->zNeutVPast = mats->oneGroupXS->zNeutV;  
  mats->oneGroupXS->rNeutVPast = mats->oneGroupXS->rNeutV;  

};
//==============================================================================

//==============================================================================
/// Run transient with multiple solves 
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

//==============================================================================
/// Run transient with multiple solves 
///
void MultiPhysicsCoupledQD::solveTransient()
{

  for (int iT = 0; iT < mesh->dts.size(); iT++)
  {
    buildLinearSystem();
    solveLinearSystem();   
    updateVarsAfterConvergence();
  }
};
//==============================================================================

//==============================================================================
/// Check for optional input parameters of relevance to this object
void MultiPhysicsCoupledQD::checkOptionalParams()
{
  string precondInput;
  
  if ((*input)["parameters"]["epsMPQD"])
  {
    epsMPQD=(*input)["parameters"]["epsMPQD"].as<double>();
  }

  if ((*input)["parameters"]["preconditionerELOT"])
  {
    precondInput=(*input)["parameters"]["preconditionerELOT"].as<string>();

    if (precondInput == "ilu")
      preconditioner = iluPreconditioner;
    else if (precondInput == "diagonal" or precondInput == "diag")
      preconditioner = diagPreconditioner;
      
  }
}
//==============================================================================
