// File: MultiGroupDNP.cpp     
// Purpose: Container for all precursor groups
// Date: April 9, 2020

#include "MultiGroupDNP.h"
#include "MultiPhysicsCoupledQD.h"
#include "SingleGroupDNP.h"

using namespace std;

//==============================================================================
/// MultiGroupDNP class object constructor
///
/// @param [in] blankType blank for this material
MultiGroupDNP::MultiGroupDNP(Materials * myMats,\
  Mesh * myMesh,\
  YAML::Node * myInput,\
  MultiPhysicsCoupledQD * myMPQD,\
  int myIndexOffset)
{
  
  // Assign pointers
  mats = myMats;
  mesh = myMesh; 
  input = myInput;
  mpqd = myMPQD;
  indexOffset = myIndexOffset;
  
  // Read input parameters
  readInput();

};
//==============================================================================

//==============================================================================
/// readInput Read in input parameters for DNPs
///
void MultiGroupDNP::readInput()
{
  // Temporary integers to store number of group unknowns
  int nGroupCoreUnknowns,nGroupRecircUnknowns,coreIndexOffset,recircIndexOffset;
  
  // Temporary vector for beta and lambda input params
  vector<double> betaInp,lambdaInp;
  Eigen::VectorXd betas,lambdas;

  // Check if DNP data are specified in input
  if ((*input)["delayed neutron precursors"]["betas"] and\
    (*input)["delayed neutron precursors"]["lambdas"])
  {
    // ToDo: allow specification of just betas or lambdas 

    betaInp = (*input)["delayed neutron precursors"]["betas"]\
            .as<vector<double>>();
    betas.setZero(betaInp.size());

    lambdaInp = (*input)["delayed neutron precursors"]["lambdas"]\
            .as<vector<double>>();
    lambdas.setZero(lambdaInp.size());

    cout << "read inputs" << endl;
    cout << "betas size: " <<  betaInp.size()<<endl;
    cout << "lambdas size: " <<  lambdaInp.size()<<endl;
    
    for (int iGroup = 0; iGroup < betaInp.size(); iGroup++)
    {
      betas(iGroup) = betaInp[iGroup];
      lambdas(iGroup) = lambdaInp[iGroup];
    } 

    beta = betas.sum(); 

  } else
  {
    // Set default DNP data.
    betas.setZero(6);
    lambdas.setZero(6);
    betas << 0.00021,0.00142,0.00128,0.00257,0.00075,0.00027;
    lambdas << 0.012375,0.03013,0.111774,0.301304,1.13607,3.01304;
    beta = betas.sum();
  }
 

  // Initialize single group DNP objects 
  nGroupCoreUnknowns = mesh->nZ*mesh->nR;
  nGroupRecircUnknowns = mesh->nZrecirc*mesh->nR;
  
  for (int iGroup = 0; iGroup < betas.size(); ++iGroup){
    coreIndexOffset = indexOffset + iGroup*nGroupCoreUnknowns;
    recircIndexOffset = iGroup*nGroupRecircUnknowns;
    shared_ptr<SingleGroupDNP> SGDNP (new SingleGroupDNP(mats,mesh,this,\
      betas(iGroup),lambdas(iGroup),coreIndexOffset,recircIndexOffset));
    DNPs.push_back(std::move(SGDNP));
  }
  
  nCoreUnknowns = DNPs.size()*nGroupCoreUnknowns;
  nRecircUnknowns = DNPs.size()*nGroupRecircUnknowns;

  // Set sizes of matrices in recirculation solve
  recircA.resize(nRecircUnknowns,nRecircUnknowns);     
  recircb.resize(nRecircUnknowns);     
  recircx.resize(nRecircUnknowns);     
};
//==============================================================================

//==============================================================================
/// Build linear system for DNP recirculation
///
void MultiGroupDNP::buildRecircLinearSystem()
{
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    DNPs[iGroup]->buildRecircLinearSystem();
  }
};
//==============================================================================

//==============================================================================
/// Build linear system for DNPs inmultiphysics coupled quasidiffusion system
///
void MultiGroupDNP::buildCoreLinearSystem()
{
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    DNPs[iGroup]->buildCoreLinearSystem();
  }
};
//==============================================================================


//==============================================================================
/// Solve linear system for multiphysics coupled quasidiffusion system
///
void MultiGroupDNP::solveRecircLinearSystem()
{
  Eigen::SparseLU<Eigen::SparseMatrix<double>,\
    Eigen::COLAMDOrdering<int> > solverLU;
  recircA.makeCompressed();
  solverLU.compute(recircA);
  recircx = solverLU.solve(recircb);
};
//==============================================================================
