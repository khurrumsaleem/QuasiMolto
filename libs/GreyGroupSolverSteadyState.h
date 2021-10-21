#ifndef GREYGROUPSOLVERSTEADYSTATE_H
#define GREYGROUPSOLVERSTEADYSTATE_H 

#include "GreyGroupSolverBase.h"

using namespace std; 

//==============================================================================
//! GreyGroupSolverSteadyState class that solves RZ neutron transport

class GreyGroupSolverSteadyState: public GreyGroupSolverBase
{
  public:

    // Constructor
    GreyGroupSolverSteadyState(GreyGroupQD * myGGQD,\
        Mesh * myMesh,\
        Materials * myMaterials,\
        YAML::Node * myInput)\
        :GreyGroupSolver(myGGQD,\
        myMesh,\ 
        myMaterials,\ 
        myInput);

};

//==============================================================================

#endif
