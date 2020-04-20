// File: SingleGroupPrecursor.cpp     
// Purpose: Contains information and builds linear system for a precursor group
// Date: April 9, 2020

#include "SingleGroupPrecursor.h"
#include "MultiGroupPrecursor.h"

using namespace std;

//==============================================================================
/// SingleGroupPrecursor class object constructor
///
/// @param [in] blankType blank for this material
SingleGroupPrecursor::SingleGroupPrecursor(MultiGroupPrecursor * myMGP,\
  double myBeta,\
  double myLambda)
{

  // Assign inputs to their member variables
  dnps = myMGP;
  beta = myBeta;
  lambda = myLambda;  
};
//==============================================================================
