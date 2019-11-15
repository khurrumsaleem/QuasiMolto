#ifndef MATERIAL_H
#define MATERIAL_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
//! Material class that holds material and geometry information

class Material
{
        public:
	int matID;
	string name;
 	Eigen::VectorXd sigT,sigF,chiP,chiD;
	Eigen::MatrixXd sigS;
	double nu;	
        // public functions
	Material(int myMatID,\
		string name,\
		Eigen::VectorXd mySigT,\
		Eigen::MatrixXd mySigS,\
		Eigen::VectorXd mySigF,\
                Eigen::VectorXd myChiP,\
                Eigen::VectorXd myChiD,\
		double myNu);
        void checkMat();
	void edit();

};

//==============================================================================

#endif
