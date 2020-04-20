// File: MultiGroupPrecursor.cpp     
// Purpose: Container for all precursor groups
// Date: April 9, 2020

#include "MultiGroupPrecursor.h"
#include "MultiPhysicsCoupledQD.h"
#include "SingleGroupPrecursor.h"

using namespace std;

//==============================================================================
/// MultiGroupPrecursor class object constructor
///
/// @param [in] blankType blank for this material
MultiGroupPrecursor::MultiGroupPrecursor(Materials * myMats,\
Mesh * myMesh,\
YAML::Node * myInput,\
MultiPhysicsCoupledQD * myMPQD)
{
  
  // Assign pointers
  mats = myMats;
  mesh = myMesh; 
  input = myInput;
  mpqd = myMPQD;
  
  // Read input parameters
  readInput();

};
//==============================================================================

//==============================================================================
/// readInput Read in input parameters for DNPs
///
void MultiGroupPrecursor::readInput()
{
  
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
  
  for (int iGroup = 0; iGroup < betas.size(); ++iGroup){
    shared_ptr<SingleGroupPrecursor> SGDNP (new SingleGroupPrecursor(this,\
      betas(iGroup),lambdas(iGroup)));
    DNPs.push_back(std::move(SGDNP));
  }
   
};
//==============================================================================
