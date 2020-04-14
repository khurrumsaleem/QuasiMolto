// File: HeatTransfer.cpp     
// Purpose: Contains information and builds linear system for heat transfer
// Date: April 9, 2020

#include "HeatTransfer.h"

using namespace std;

//==============================================================================
/// HeatTransfer class object constructor
///
/// @param [in] myMaterials materials object for the simulation
/// @param [in] myMesh mesh object for the simulation
/// @param [in] myInput input object for the simulation
/// @param [in] myQD multiphysics coupled QD object for the simulation
HeatTransfer::HeatTransfer(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput,\
  MultiPhysicsCoupledQD * myQD)
{
  // Assign inputs to their member variables
  mats = myMaterials;
  mesh = myMesh;
  input = myInput;
  qd = myQD;

  // Initialize size of matrices
  temp.setZero(mesh->nZ,mesh->nR);
  flux.setZero(mesh->nZ+1,mesh->nR);
  dirac.setZero(mesh->nZ+1,mesh->nR);

  // Check for optional inputs 
  if ((*input)["parameters"]["wallTemp"]){
    wallT=(*input)["parameters"]["wallTemp"].as<double>();
  } 
  if ((*input)["parameters"]["inletTemp"]){
    inletT=(*input)["parameters"]["inletTemp"].as<double>();
  } 

  cout << "Initialized HeatTransfer object." << endl;
  cout << "wallTemp: " << wallT << endl;
  cout << "inletTemp: " << inletT << endl;
  cout << getIndex(4,4) << endl; 
};
//==============================================================================

//==============================================================================
/// Calculate energy diracs
///
void HeatTransfer::calcDiracs()
{
  bool posVelocity = true;
  double TupwindInterface,Tinterface,theta,phi;
  Eigen::Vector2d fluxLimiterArg1,fluxLimiterArg2;
  Eigen::Vector3d fluxLimiterArg3;

  if (posVelocity) 
  {
    for (int iR = 0; iR < dirac.cols()-1; iR++)
    {
      // Handle iZ=1 case

      // Handle iZ=2 case
      
      // Handle all other cases
      for (int iZ = 2; iZ < dirac.rows()-1; iZ++)
      {

        TupwindInterface = temp(iZ-1,iR) - temp(iZ-2,iR);
        Tinterface = temp(iZ,iR) - temp(iZ-1,iR);
        theta = TupwindInterface/Tinterface;
        fluxLimiterArg1 << 1,2*theta; 
        fluxLimiterArg2 << 2,theta; 
        fluxLimiterArg3 << 0,\
          fluxLimiterArg1.minCoeff(),\
          fluxLimiterArg2.minCoeff();
        phi = fluxLimiterArg3.maxCoeff();
        dirac(iZ,iR) = phi*Tinterface; 
        
      }
    }
  } else 
  {

  }
  
};
//==============================================================================


//==============================================================================
/// Calculate energy flux 
///
void HeatTransfer::calcFluxes()
{
  bool posVelocity = true;

  if (posVelocity) {
    
  }

  
};
//==============================================================================


//==============================================================================
/// Map 2D coordinates to index of temperature in the 1D solution vector
///
/// @param [in] iZ axial location
/// @param [in] iR radial location
/// @param [out] index the index for temperature in the 1D solution vector
int HeatTransfer::getIndex(int iZ, int iR)
{

  int index,nR = mesh->drsCorner.size();

  index = indexOffset + iR + nR*iZ;

  return index;
  
};
//==============================================================================

//==============================================================================
/// Map 2D coordinates to index of temperature in the 1D solution vector
///
/// @param [in] iZ axial location
/// @param [in] iR radial location
/// @param [out] index the index for temperature in the 1D solution vector
void HeatTransfer::getTemp()
{

  for (int iR = 0; iR < mesh-> drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh-> dzsCorner.size(); iZ++)
    { 
    
      temp(iZ,iR) = qd->x(getIndex(iZ,iR));   
    
    }
  }
  
};
//==============================================================================
