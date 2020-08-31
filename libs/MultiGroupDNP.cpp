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
  double tempVal;  


  // Temporary vector for beta and lambda input params
  vector<double> betaInp,lambdaInp;
  Eigen::VectorXd lambdas;
  Eigen::MatrixXd betas;

  // Set size of beta vector
  beta.setZero(mats->nGroups);
  
  // set size of cumulative DNP decay source 
  dnpSource.setZero(mesh->nZ,mesh->nR);

  // Check if DNP data are specified in input
  if ((*input)["delayed neutron precursors"]["betas"] and\
      (*input)["delayed neutron precursors"]["lambdas"])
  {
    // ToDo: allow specification of just betas or lambdas 

    lambdaInp = (*input)["delayed neutron precursors"]["lambdas"]\
                .as<vector<double>>();
    lambdas.setZero(lambdaInp.size());

    betaInp = (*input)["delayed neutron precursors"]["betas"]\
              .as<vector<double>>();
    betas.setZero(mats->nGroups,lambdaInp.size());

    for (int iDNPGroup = 0; iDNPGroup < lambdaInp.size(); iDNPGroup++)
    {
      lambdas(iDNPGroup) = lambdaInp[iDNPGroup];
    } 
    // If the size of the beta input matches that of the lambda input, then
    // assume the same set of betas is used in each neutron energy group
    if (lambdaInp.size() == betaInp.size())
    {
      for (int iDNPGroup = 0; iDNPGroup < betaInp.size(); iDNPGroup++)
      {
        for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
          betas(iEnergyGroup,iDNPGroup) = betaInp[iDNPGroup];
      }
    }
    // Otherwise, assume betas are specified for each neutron energy group and 
    // cycle through the betaInp appropriately 
    else
    {
      for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
      {
        for (int iDNPGroup = 0; iDNPGroup < lambdaInp.size(); iDNPGroup++)
        {
          betas(iEnergyGroup,iDNPGroup) 
            = betaInp[iEnergyGroup*lambdaInp.size() + iDNPGroup];
        }
      }
    } 
  } 
  // If betas and lambdas are not specified, assume the defaults defined below
  else
  {
    // Set default DNP data.
    betas.setZero(mats->nGroups,6);
    lambdas.setZero(6);
    for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
    {
      betas.row(iEnergyGroup) << 0.00021,0.00142,0.00128,0.00257,0.00075,\
        0.00027;
    }
    lambdas << 0.012375,0.03013,0.111774,0.301304,1.13607,3.01304;

  }
  // Calculate total DNP fraction in each neutron energy group
  for (int iEnergyGroup = 0; iEnergyGroup < mats->nGroups; iEnergyGroup++)
  {     
    beta(iEnergyGroup) = betas.row(iEnergyGroup).sum();
  }

  // Initialize single group DNP objects 
  nGroupCoreUnknowns = mesh->nZ*mesh->nR;
  nGroupRecircUnknowns = mesh->nZrecirc*mesh->nR;

  for (int iDNPGroup = 0; iDNPGroup < lambdas.size(); ++iDNPGroup){
    coreIndexOffset = indexOffset + iDNPGroup*nGroupCoreUnknowns;
    recircIndexOffset = iDNPGroup*nGroupRecircUnknowns;
    shared_ptr<SingleGroupDNP> SGDNP (new SingleGroupDNP(mats,mesh,this,\
          betas.col(iDNPGroup),lambdas(iDNPGroup),coreIndexOffset,\
          recircIndexOffset,iDNPGroup));
    DNPs.push_back(std::move(SGDNP));
  }

  nCoreUnknowns = DNPs.size()*nGroupCoreUnknowns;
  nRecircUnknowns = DNPs.size()*nGroupRecircUnknowns;

  // Set sizes of matrices in recirculation solve
  recircA.resize(nRecircUnknowns,nRecircUnknowns);     
  recircb.resize(nRecircUnknowns);     
  recircx.resize(nRecircUnknowns);    

  // Initialize DNP data in collapsed cross section object
  mats->oneGroupXS->groupDNPFluxCoeff.resize(lambdas.size()); 

  for (int iDNPGroup = 0; iDNPGroup < lambdas.size(); iDNPGroup ++)
  {
    mats->oneGroupXS->groupDNPFluxCoeff[iDNPGroup].setZero(mesh->nZ,mesh->nR);

    for (int iR = 0; iR < mesh->nR; iR++)
    {
      for (int iZ = 0; iZ < mesh->nZ; iZ++)
      {

        tempVal = DNPs[iDNPGroup]->beta(0)*mats->nu(iZ,iR,0)\
                  *mats->oneGroupXS->sigF(iZ,iR);  
        mats->oneGroupXS->groupDNPFluxCoeff[iDNPGroup](iZ,iR) = tempVal;

      }
    } 
  }

  for (int iR = 0; iR < mesh->nR; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ; iZ++)
    {

      tempVal = (1-beta(0))*mats->nu(iZ,iR,0)\
                *mats->oneGroupXS->sigF(iZ,iR);  
      mats->oneGroupXS->qdFluxCoeff(iZ,iR) = tempVal;

    }
  } 
};
//==============================================================================

//==============================================================================
/// Build linear system for DNP recirculation
///
void MultiGroupDNP::buildRecircLinearSystem()
{
  recircA.setZero();
  recircb.setZero();
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    DNPs[iGroup]->buildRecircLinearSystem();
  }
};
//==============================================================================

//==============================================================================
/// Build linear system for steady state DNP recirculation
///
void MultiGroupDNP::buildSteadyStateRecircLinearSystem()
{
  recircA.setZero();
  recircb.setZero();
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    DNPs[iGroup]->buildSteadyStateRecircLinearSystem();
  }
};
//==============================================================================

//==============================================================================
/// Build linear system for transient DNPs in multiphysics coupled 
/// quasidiffusion system
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
/// Build linear system for steady state DNPs in multiphysics coupled 
/// quasidiffusion system
///
void MultiGroupDNP::buildSteadyStateCoreLinearSystem()
{
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    DNPs[iGroup]->buildSteadyStateCoreLinearSystem();
  }
};
//==============================================================================


//==============================================================================
/// Extract core DNP concentrations into 2D matrices in each group 
///
void MultiGroupDNP::getCoreDNPConc()
{
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    DNPs[iGroup]->getCoreConc();
  }

  getCumulativeDNPDecaySource();
};
//==============================================================================

//==============================================================================
/// Update matrix holding cumulative DNP decay source 
///
void MultiGroupDNP::getCumulativeDNPDecaySource()
{
  double groupLambda,localGroupConc;
  int groupOffset; 
  dnpSource.setZero();
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    groupLambda = DNPs[iGroup]->lambda;
    groupOffset = DNPs[iGroup]->coreIndexOffset;
    for (int iR = 0; iR < mesh->nR; iR++)
    {
      for (int iZ = 0; iZ < mesh->nZ; iZ++)
      {
        localGroupConc = mpqd->x(DNPs[iGroup]->getIndex(iZ,iR,groupOffset)); 
        dnpSource(iZ,iR) += groupLambda*localGroupConc;  
      }
    }
    //dnpSource += DNPs[iGroup]->lambda*DNPs[iGroup]->dnpConc;
  }
};
//==============================================================================

//==============================================================================
/// Extract recirc DNP concentrations into 2D matrices in each group 
///
void MultiGroupDNP::getRecircDNPConc()
{
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    DNPs[iGroup]->getRecircConc();
  }
};
//==============================================================================

//==============================================================================
/// Extract core DNP concentrations into 2D matrices in each group 
///
void MultiGroupDNP::setCoreDNPConc()
{
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    DNPs[iGroup]->setCoreConc();
  }
};
//==============================================================================

//==============================================================================
/// Extract core DNP concentrations into 2D matrices in each group 
///
void MultiGroupDNP::setRecircDNPConc()
{
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    DNPs[iGroup]->setRecircConc();
  }
};
//==============================================================================

//==============================================================================
/// Print core DNP concentrations in 2D matrices 
///
void MultiGroupDNP::printCoreDNPConc()
{
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    cout << "Group: " << iGroup << endl;
    cout << DNPs[iGroup]->dnpConc << endl;
  }
};
//==============================================================================

//==============================================================================
/// Print core DNP concentrations in 2D matrices 
///
void MultiGroupDNP::printRecircDNPConc()
{
  for (int iGroup = 0; iGroup < DNPs.size(); ++iGroup)
  {
    cout << "Group: " << iGroup << endl;
    cout << DNPs[iGroup]->recircConc << endl;
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
