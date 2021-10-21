#ifndef GREYGROUPSOLVERTRANSIENT_H
#define GREYGROUPSOLVERTRANSIENT_H 

#include "GreyGroupSolverBase.h"

using namespace std; 

//==============================================================================
//! GreyGroupSolver class that solves RZ neutron transport

class GreyGroupSolverTransient: public GreyGroupSolverBase
{
  public:

    // Constructor
    GreyGroupSolverTransient(GreyGroupQD * myGGQD,\
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
