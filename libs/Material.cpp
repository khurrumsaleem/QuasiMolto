// File: Material.cpp
// Purpose: define a class that holds information for each material
// Date: October 28, 2019

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Material.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
//! Material class object constructor

Material::Material(int myMatID,\
	string myName,\
	Eigen::VectorXd mySigT,\
	Eigen::MatrixXd mySigS,\
	Eigen::VectorXd mySigF,\
	double myNu)
{
	matID = myMatID;
	name = myName;
	sigT = mySigT;
	sigS = mySigS;
	sigF = mySigF;
	nu = myNu;
};

//==============================================================================

//==============================================================================
//! edit Print material information

void Material::edit()
{
int spacing = 5;
int spacing2 = 5;
cout << endl;
cout << "==========================================="<< endl;
cout << "MATERIAL: " << setw(spacing2) << name << "\n";
cout << "-------------------------------------------"<< endl;
cout << "matID:"<< matID << endl;
cout << endl;
cout << "sigT:"<<endl;
cout << sigT.transpose() << endl;
cout << endl;
cout << "sigS:"<< endl;
cout << sigS << endl;
cout << endl;
cout << "sigF:"<< endl; 
cout << sigF.transpose() << endl;
cout << endl;
cout << "nu:  "<< nu << endl;
cout << "==========================================="<< endl;
};

//==============================================================================


