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
  ggqd->GGSolver->assignPointers(&A,&x,&xPast,&b);

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
void MultiPhysicsCoupledQD::fluxSource(int iZ,int iR,int iEq,double coeff)
{

  int iCF = 0; // index of cell-average flux value in index vector  
  vector<int> indices = ggqd->GGSolver->getIndices(iR,iZ);

  A.coeffRef(iEq,indices[0]) += coeff; 

};
//==============================================================================

//==============================================================================
/// Include a dnp source in the linear system
///
/// @param [in] iZ axial location 
/// @param [in] iR radial location
/// @param [in] iEq equation index
/// @param [in] coeff coefficient of dnp source
void MultiPhysicsCoupledQD::dnpSource(int iZ,int iR,int iEq,double coeff)
{
  int index,indexOffset;
  double groupLambda;

  for (int iGroup = 0; iGroup < mgdnp->DNPs.size(); ++iGroup)
  {
    indexOffset = mgdnp->DNPs[iGroup]->coreIndexOffset;
    index = mgdnp->DNPs[iGroup]->getIndex(iZ,iR,indexOffset);
    groupLambda = mgdnp->DNPs[iGroup]->lambda;

    A.coeffRef(iEq,index) += coeff*groupLambda;
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
  
 // Eigen::SparseLU<Eigen::SparseMatrix<double>,\
    Eigen::COLAMDOrdering<int> > solverLU;
  Eigen::SuperLU<Eigen::SparseMatrix<double>> solverLU;
  A.makeCompressed();
  solverLU.compute(A);
  x = solverLU.solve(b);

  mgdnp->solveRecircLinearSystem();

};
//==============================================================================

//==============================================================================
/// Solve linear system for multiphysics coupled quasidiffusion system with a 
/// solver that can us multiple processors
///
void MultiPhysicsCoupledQD::solveLinearSystemParallel()
{

  Eigen::MatrixXd A_dense;
  A_dense = A;
  x = A_dense.partialPivLu().solve(b);

  mgdnp->solveRecircLinearSystem();

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
  if ((*input)["parameters"]["epsMPQD"])
  {
    epsMPQD=(*input)["parameters"]["epsMPQD"].as<double>();
  }
}
//==============================================================================
