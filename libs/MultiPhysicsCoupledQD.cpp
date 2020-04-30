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
  int index,groupBeta,groupLambda,indexOffset;

  for (int iGroup = 0; iGroup < mgdnp->DNPs.size(); ++iGroup)
  {
    indexOffset = mgdnp->DNPs[iGroup]->coreIndexOffset;
    index = mgdnp->DNPs[iGroup]->getIndex(iZ,iR,indexOffset);
    groupBeta = mgdnp->DNPs[iGroup]->beta;
    groupLambda = mgdnp->DNPs[iGroup]->lambda;

    A.coeffRef(iEq,index) += coeff*groupLambda*groupBeta;
  }

};
//==============================================================================

//==============================================================================
/// Build linear system for multiphysics coupled quasidiffusion system
///
void MultiPhysicsCoupledQD::buildLinearSystem()
{
  // Reset linear system
  A.setZero();
  x.setZero();
  b.setZero();
  cout << "Reset linear system" << endl; 

  // Build QD system
  ggqd->buildLinearSystem();
  cout << "Built GGQD" << endl; 
  
  // Build heat transfer system
  heat->buildLinearSystem();
  cout << "Built heat" << endl; 

  // Build delayed neutron precursor balance system
  mgdnp->buildCoreLinearSystem();  
  cout << "Built DNP" << endl; 

};
//==============================================================================

//==============================================================================
/// Solve linear system for multiphysics coupled quasidiffusion system
///
void MultiPhysicsCoupledQD::solveLinearSystem()
{
  Eigen::SparseLU<Eigen::SparseMatrix<double>,\
    Eigen::COLAMDOrdering<int> > solverLU;
  A.makeCompressed();
  solverLU.compute(A);
  x = solverLU.solve(b);

  ggqd->GGSolver->getFlux();
  cout << "flux" << endl; 
  cout << ggqd->sFlux << endl; 
  
  heat->getTemp();
  cout << "temp" << endl; 
  cout << heat->temp << endl; 
  
  mgdnp->getCoreDNPConc();
  cout << "DNP conc" << endl; 
  mgdnp->printCoreDNPConc();
};
//==============================================================================
