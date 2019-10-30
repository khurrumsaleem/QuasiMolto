#ifndef MATERIALS_H
#define MATERIALS_H

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
	umat matMap; //resized later
        // public functions
        Materials(Mesh * myMesh,YAML::Node * myInput);
	void readMats();
	void readGeom();
        void setMatRegion(int myIndex,double rIn,double rOut,\
		double zUp,double zLow);

        private:
        // private functions
        YAML::Node * input;
        Mesh * mesh;
};

//==============================================================================

#endif
