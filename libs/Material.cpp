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
	Eigen::VectorXd myChiP,\
	Eigen::VectorXd myChiD,\
	double myNu)
{
	matID = myMatID;
	name = myName;
	sigT = mySigT;
	sigS = mySigS;
	sigF = mySigF;
        chiP = myChiP;
        chiD = myChiD;
	nu = myNu;
};

//==============================================================================

//==============================================================================
//! checkMat Make sure that the input data makes sense

void Material::checkMat()
{
  double totalSigS=0.0,totalChiD=0.0,totalChiP=0.0;
 
  for (int iGroup = 0; iGroup < sigS.rows(); ++iGroup){ 

    totalSigS=0;
  
    // take sum of group transfer cross sections
    for (int iSig = 0; iSig < sigS.cols(); ++iSig){
    
      totalSigS = totalSigS + sigS(iGroup,iSig);

    }

    // check to see if the sum of the fission and scattering cross sections
    // exceed the total cross section
    if (totalSigS+sigF(iGroup)>sigT(iGroup)){
      
      cout << "***WARNING***: Material '" << name <<"', group " << iGroup <<endl;
      cout<< "Sum of scattering and fission cross sections ";
      cout << "exceeds total cross section." << endl;
      cout << endl;
      cout << "sigT:"<<endl;
      cout << sigT(iGroup) << endl;
      cout << endl;
      cout << "sigS:"<< endl;
      cout << sigS.row(iGroup) << endl;
      cout << endl;
      cout << "sigF:"<< endl; 
      cout << sigF(iGroup) << endl;
      cout << endl;

    }
    totalChiP = totalChiP + chiP(iGroup);
    totalChiD = totalChiD + chiD(iGroup);
  }

  // check to see if chiP does not sum to one
  if (totalChiP != 1.0){
      
    cout << "***WARNING***: Material '" << name <<"'"<< endl;
    cout << "ChiP does not sum to 1.0." << endl;
    cout << endl;

  }

  // check to see if chiD does not sum to one
  if (totalChiD != 1.0){
      
    cout << "***WARNING***: Material '" << name <<"'"<<endl;
    cout<< "ChiD does not sum to 1.0." << endl;
    cout << endl;

  }
};

//==============================================================================


//==============================================================================
//! edit Print material information

void Material::edit()
{
int spacing = 5;
int spacing2 = 5;
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
cout << endl;
cout << "chiP:"<< endl;
cout << chiP.transpose() << endl;
cout << endl;
cout << "chiD:"<< endl;
cout << chiD.transpose() << endl;
cout << endl;
cout << "==========================================="<< endl;
cout << endl;
};

//==============================================================================


