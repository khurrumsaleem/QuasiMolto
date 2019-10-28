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
        StartingAngle(Mesh::Mesh myMesh,YAML::Node myInput);

        private:
        // private functions
        YAML::Node * input;
        Mesh::Mesh * mesh;
};

//==============================================================================

//==============================================================================
//! StartingAngle object constructor

StartingAngle::StartingAngle(Mesh::Mesh myMesh,\
                             YAML::Node myInput)
{
	// Point to variables for mesh and input file
	mesh = &myMesh;
	input = &myInput;
};

//==============================================================================
