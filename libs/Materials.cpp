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
	umat matMap;
        // public functions
        Materials(Mesh * myMesh,YAML::Node * myInput);
	void readInput();
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
	matMap.set_size(mesh->zCent.size(),mesh->rCent.size());
	matMap.fill(0);
	readInput();
};

//==============================================================================

//==============================================================================
//! readInput Parse input and define material and geometry parameters

void Materials::readInput()
{
	YAML::Node geometry = (*input)["geometry"];
	string regType;
        
	cout << "size of geometry" << geometry.size() << endl;
	for (YAML::const_iterator it=geometry.begin(); it!=geometry.end(); ++it){
		regType=it->first.as<string>();
		if (regType=="background"){
			cout << "blanke region in type!" << endl;
		} else if (regType=="region"){
			cout << "fill specific region"<< endl;
			cout << it->second["material"].as<string>() << endl;
		}
	}
};

//==============================================================================


