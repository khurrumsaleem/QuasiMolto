// File: SingleGroupDNP.cpp     
// Purpose: Contains information and builds linear system for a precursor group
// Date: April 9, 2020

#include "SingleGroupDNP.h"
#include "MultiGroupDNP.h"

using namespace std;

//==============================================================================
/// SingleGroupDNP class object constructor
///
/// @param [in] blankType blank for this material
SingleGroupDNP::SingleGroupDNP(MultiGroupDNP * myMGDNP,\
  double myBeta,\
  double myLambda)
{

  // Assign inputs to their member variables
  dnps = myMGDNP;
  beta = myBeta;
  lambda = myLambda;  
};
//==============================================================================
