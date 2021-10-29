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
        :GreyGroupSolverBase(myGGQD,\
        myMesh,\
        myMaterials,\
        myInput){};

    // functions to enforce governing equations
    int assertZerothMoment(int iR,int iZ,int iEq);

    // functions to enforce coefficients for facial currents
    int southCurrent(double coeff,int iR,int iZ,int iEq);
    int northCurrent(double coeff,int iR,int iZ,int iEq);
    int westCurrent(double coeff,int iR,int iZ,int iEq);
    int eastCurrent(double coeff,int iR,int iZ,int iEq);

    // functions to enforce coefficients for calculation of facial currents
    int calcSouthCurrent(int iR,int iZ,int iEq);
    int calcNorthCurrent(int iR,int iZ,int iEq);
    int calcWestCurrent(int iR,int iZ,int iEq);
    int calcEastCurrent(int iR,int iZ,int iEq);

};

//==============================================================================

#endif
