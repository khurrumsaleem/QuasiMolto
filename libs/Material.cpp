// File: Material.cpp
// Purpose: define a class that holds information for each material
// Date: October 28, 2019

#include "Material.h"

using namespace std; 

//==============================================================================
/// Material class object constructor
///
/// @param [in] myMatID Material index for this material
/// @param [in] myName Name of this material
/// @param [in] mySigT Total cross sections for this material
/// @param [in] mySigS Scattering cross sections for this material
/// @param [in] mySigF Fission cross sections for this material
/// @param [in] myChiP Prompt fission probabilities in this material 
/// @param [in] myChiD Delayed fission probabilities in this material 
/// @param [in] nu Average number of neutrons born per fission event in this 
/// material
Material::Material(int myMatID,\
  string myName,\
  vector<Eigen::MatrixXd> mySigT,\
  vector<Eigen::MatrixXd> mySigS,\
  vector<Eigen::MatrixXd> mySigF,\
  vector<Eigen::MatrixXd> myNu,\
  vector<Eigen::MatrixXd> myNeutV,\
  Eigen::MatrixXd myFlowVelocity,\
  Eigen::VectorXd myChiP,\
  Eigen::VectorXd myChiD,\
  double myDensity,\
  double myGamma,\
  double myK,\
  double mycP,\
  double myOmega,\
  bool myStationary)
{
  // Assign inputs to their member variables
  matID = myMatID;
  name = myName;
  sigT = mySigT;
  sigS = mySigS;
  sigF = mySigF;
  chiP = myChiP;
  chiD = myChiD;
  nu = myNu;
  neutV = myNeutV;
  density = myDensity;
  gamma = myGamma;
  k = myK;
  cP = mycP;
  omega = myOmega;
  flowVelocity = myFlowVelocity;
  stationary = myStationary;

};
//==============================================================================

//==============================================================================
/// Return total cross section 
/// 
/// @param [in] eIdx energy index for cross section 
/// @param [in] temp temperature to evaluate XS at 
/// @param [out] interpolated XS 
double Material::getSigT(int eIdx,double temp)
{
  return interpolateParameter(sigT[eIdx],temp); 
};

//==============================================================================

//==============================================================================
/// Return fission cross section 
/// 
/// @param [in] eIdx energy index for cross section 
/// @param [in] temp temperature to evaluate XS at 
/// @param [out] interpolated XS 
double Material::getSigF(int eIdx,double temp)
{
  return interpolateParameter(sigF[eIdx],temp); 
};

//==============================================================================

//==============================================================================
/// Return fission cross section 
/// 
/// @param [in] eIdxPrim from energy group
/// @param [in] eIdx to energy group
/// @param [in] temp temperature to evaluate XS at 
/// @param [out] interpolated XS 
double Material::getSigS(int eIdxPrime,int eIdx,double temp)
{
  int size = sigT.size();
  return interpolateParameter(sigS[size*eIdxPrime+eIdx],temp); 
};

//==============================================================================

//==============================================================================
/// Return nu 
/// 
/// @param [in] eIdx energy index for nu
/// @param [in] temp temperature to evaluate nu at 
/// @param [out] interpolated nu 
double Material::getNu(int eIdx,double temp)
{
  return interpolateParameter(nu[eIdx],temp); 
};

//=============================================================================

//=============================================================================
/// Return neutron velocity 
/// 
/// @param [in] eIdx energy index for neutron velocity
/// @param [in] temp temperature to evaluate neutron velocity at 
/// @param [out] interpolated neutron velocity 
double Material::getNeutV(int eIdx,double temp)
{
  return interpolateParameter(neutV[eIdx],temp); 
};

//==============================================================================

//=============================================================================
/// Return flow velocity 
/// 
/// @param [in] time to evaluate neutron velocity at 
/// @param [out] interpolated flow velocity 
double Material::getFlowVelocity(double time)
{
  return interpolateParameter(flowVelocity,time); 
};

//==============================================================================

//==============================================================================
/// Interpolates parameters 
/// 
/// @param [in] param list of temperature dependend params 
/// @param [in] temp temperature to evaluate XS at 
/// @param [out] interpolated parameter 
double Material::interpolateParameter(Eigen::MatrixXd param,double temp)
{
  
  double lowerIndex,upperIndex;
  double lowerTemp,upperTemp;
  double lowerParam,upperParam;
  double slope,interpParam;
  double eps = 1E-6;
  int lastIndex = param.rows()-1;

  // Check if temperature is lower than specified in table or if 
  // there's only one entry in the temperature table
  if (temp - eps < param(0,0) or param.rows() == 1)
  {
    // Just use lower bound of parameter
    interpParam = param(0,1); 
  }
  // Check if temperature is greater than specified in table
  else if (temp + eps > param(lastIndex,0))
  {
    // Just use upper bound of parameter
    interpParam = param(lastIndex,1); 
  }
  // Otherwise, interpolate
  else
  {
    // Check if temperature is greater than specified in table 
    for (int iRow = 0; iRow < param.rows(); iRow++)
    {
      if (temp + eps < param(iRow + 1,0))
      {
        lowerTemp = param(iRow, 0); 
        lowerParam = param(iRow, 1); 
        upperTemp = param(iRow + 1, 0); 
        upperParam = param(iRow + 1, 1); 
        break; 
      }
    }
    // Linear interpolation
    slope = (upperParam - lowerParam)/(upperTemp-lowerTemp);
    interpParam = lowerParam + slope * (temp - lowerTemp); 
  }

  return interpParam;
};
//==============================================================================

//==============================================================================
/// Make sure that the input data for this material makes sense

void Material::checkMat()
{
  double totalSigS=0.0,totalChiD=0.0,totalChiP=0.0;

//  // Loop over groups 
//  for (int iGroup = 0; iGroup < sigS.rows(); ++iGroup){ 
//
//    totalSigS=0;
//  
//    // Take sum of group transfer cross sections
//    for (int iSig = 0; iSig < sigS.cols(); ++iSig){
//    
//      totalSigS = totalSigS + sigS(iGroup,iSig);
//
//    }
//
//    // Check to see if the sum of the fission and scattering cross sections
//    // exceed the total cross section
//    if (totalSigS+sigF(iGroup)>sigT(iGroup)){
//     
//      // Print a warning 
//      cout << "***WARNING***: Material '" << name <<"', group " << iGroup <<endl;
//      cout<< "Sum of scattering and fission cross sections ";
//      cout << "exceeds total cross section." << endl;
//      cout << endl;
//      cout << "sigT:"<<endl;
//      //cout << sigT[iGroup] << endl;
//      cout << endl;
//      cout << "sigS:"<< endl;
//      //cout << sigS[iGroup].row(iGroup) << endl;
//      cout << endl;
//      cout << "sigF:"<< endl; 
//      //cout << sigF(iGroup) << endl;
//      cout << endl;
//
//    }
//    // Accumulate chiPs and chiDs
//    totalChiP = totalChiP + chiP(iGroup);
//    totalChiD = totalChiD + chiD(iGroup);
//  }
//
//  // Check to see if chiP does not sum to one
//  if (totalChiP != 1.0){
//      
//    cout << "***WARNING***: Material '" << name <<"'"<< endl;
//    cout << "ChiP does not sum to 1.0." << endl;
//    cout << endl;
//
//  }
//
//  // Check to see if chiD does not sum to one
//  if (totalChiD != 1.0){
//      
//    cout << "***WARNING***: Material '" << name <<"'"<<endl;
//    cout<< "ChiD does not sum to 1.0." << endl;
//    cout << endl;
//
//  }
};

//==============================================================================


//==============================================================================
/// Print material information

void Material::edit()
{
  int spacing = 5,spacing2 = 5;
  cout << "==========================================="<< endl;
  cout << "MATERIAL: " << setw(spacing2) << name << "\n";
  cout << "-------------------------------------------"<< endl;
  cout << "matID:"<< matID << endl;
  cout << endl;
  cout << "sigT:"<<endl;
//  cout << sigT.transpose() << endl;
  cout << endl;
  cout << "sigS:"<< endl;
//  cout << sigS << endl;
  cout << endl;
  cout << "sigF:"<< endl; 
//  cout << sigF.transpose() << endl;
  cout << endl;
//  cout << "nu:  "<< nu << endl;
  cout << endl;
  cout << "chiP:"<< endl;
  cout << chiP.transpose() << endl;
  cout << endl;
  cout << "chiD:"<< endl;
  cout << chiD.transpose() << endl;
  cout << endl;
  cout << "density: "<< density << endl;
  cout << endl;
  cout << "k: "<< k << endl;
  cout << endl;
  cout << "cP: "<< cP << endl;
  cout << endl;
  cout << "gamma: "<< gamma << endl;
  cout << endl;
  cout << "omega: "<< omega << endl;
  cout << endl;
  cout << "stationary: "<< stationary << endl;
  cout << endl;
  cout << "==========================================="<< endl;
  cout << endl;
};

//==============================================================================


