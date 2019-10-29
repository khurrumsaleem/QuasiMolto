// File: StartingAngle.cpp
// Purpose: calculate initial half angle fluxes needed for 
// approximation of the angular redistribution term of the
// neutron transport equation in RZ geometry
// Date: October 28, 2019

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "Mesh.h"

using namespace std; 
using namespace arma;

//==============================================================================
//! StartingAngle class that solves RZ neutron transport at the starting angles

class StartingAngle
{
        public:
        // public functions
        StartingAngle(Mesh myMesh,YAML::Node myInput);
        void calcStartingAngle();

        private:
        // private functions
        YAML::Node * input;
        Mesh * mesh;
};

//==============================================================================

//==============================================================================
//! StartingAngle object constructor

StartingAngle::StartingAngle(Mesh myMesh,\
                             YAML::Node myInput)
{
	// Point to variables for mesh and input file
	mesh = &myMesh;
	input = &myInput;
};

//==============================================================================

//==============================================================================
//! calcStartingAngle Calculate the starting half angle angular flux

// Solves using simple corner balance on a simplified neutron transport equation 
// in RZ geometry. The simplification comes about by substituting quadrature 
// values that correspond to the starting half angle, and results in the 
// elimination of the angular redistribution term. 

void StartingAngle::calcStartingAngle()
{
        // index xi value is stored in in quadLevel
	const int xiIndex = 0;
        // temporary variable used for looping though quad set
	double xi;
        	

	for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi){
		// get value for xi at this quadrature level
		xi = mesh->quadrature[iXi].quad[0][xiIndex];
	}
	

};

//==============================================================================
