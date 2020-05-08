#ifndef MATERIAL_H
#define MATERIAL_H

#include "Mesh.h"

using namespace std; 
using namespace arma;

//==============================================================================
//! Material class that holds material and geometry information

class Material
{
        public:
	int matID;
	string name;
 	Eigen::VectorXd sigT,sigF,chiP,chiD,nu;
	Eigen::MatrixXd sigS;
        double density; // density of the material
        double gamma; // fraction of energy deposited from gamma rays
        double k; // thermal conductivity
        double cP; // specific heat
        double omega; // energy release per fission	
        bool stationary = true; // is the material stationary
        
        // public functions
	Material(int myMatID,\
		string name,\
		Eigen::VectorXd mySigT,\
		Eigen::MatrixXd mySigS,\
		Eigen::VectorXd mySigF,\
                Eigen::VectorXd myChiP,\
                Eigen::VectorXd myChiD,\
		Eigen::VectorXd myNu,\
                double myDensity,\
                double myGamma,\
                double myK,\
                double mycP,\
                double myOmega,\
                bool myStationary);
        void checkMat();
	void edit();

};

//==============================================================================

#endif
