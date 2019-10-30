// File: Materials.cpp
// Purpose: contain material characteristics and geometry for a simulation
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
//! Material class that holds material and geometry information

class Materials
{
        public:
        // public functions
        Material(Mesh * myMesh,YAML::Node * myInput);

        private:
        // private functions
        YAML::Node * input;
        Mesh * mesh;
};

//==============================================================================

//==============================================================================
//! Material class object constructor

Materials::Materials(Mesh * myMesh,\
                             YAML::Node * myInput)
{
	// Point to variables for mesh and input file
	mesh = myMesh;
	input = myInput;
};

//==============================================================================

