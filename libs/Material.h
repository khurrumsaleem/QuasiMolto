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
 	Eigen::VectorXd sigT,sigF,chiP,chiD;
	Eigen::MatrixXd sigS;
	double nu; // average number of neutrons released per fission 
        double density; // density of the material
        double gamma; // fraction of energy deposited from gamma rays
        double k; // thermal conductivity
        double cP; // specific heat
        double omega; // energy release per fission	
        
        // public functions
	Material(int myMatID,\
		string name,\
		Eigen::VectorXd mySigT,\
		Eigen::MatrixXd mySigS,\
		Eigen::VectorXd mySigF,\
                Eigen::VectorXd myChiP,\
                Eigen::VectorXd myChiD,\
		double myNu,\
                double myDensity,\
                double myGamma,\
                double myK,\
                double mycP,\
                double myOmega);
        void checkMat();
	void edit();

};

//==============================================================================

#endif
