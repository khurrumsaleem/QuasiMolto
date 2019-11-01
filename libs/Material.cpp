// File: Material.cpp
// Purpose: define a class that holds information for each material
// Date: October 28, 2019

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
	string name;
 	rowvec sigT,sigF;
	mat sigS;
	double nu;	
        // public functions
	Material(int myMatID,\
		string name,\
		rowvec mySigT,\
		mat mySigS,\
		rowvec mySigF,\
		double myNu);
	void edit();

};

//==============================================================================

//==============================================================================
//! Material class object constructor

Material::Material(int myMatID,\
	string myName,\
	rowvec mySigT,\
	mat mySigS,\
	rowvec mySigF,\
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
//! Material class object constructor

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
cout << sigT << endl;
cout << "sigS:"<< endl;
cout << sigS << endl;
cout << "sigF:"<< endl; 
cout << sigF << endl;
cout << "nu:  "<< nu << endl;
cout << "==========================================="<< endl;
};

//==============================================================================


