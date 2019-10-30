#ifndef MATERIALS_H
#define MATERIALS_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"

using namespace std; 
using namespace arma;

//==============================================================================
//! Material class that holds material and geometry information

class Material
{
        public:
	int matID;
 	rowvec sigT;
	rowvec sigS;
	rowvec sigF;
	double nu;	
        // public functions
	Material(int myMatID,\
		rowvec mySigT,\
		rowvec mySigS,\
		rowvec mySigF,\
		double myNu)

};

//==============================================================================

#endif
