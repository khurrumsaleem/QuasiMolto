// File: SimpleCornerBalance.cpp
// Purpose: Solve RZ neutron transport equation using simple corner balance 
// Date: October 28, 2019

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Materials.h"
#include "Material.h"
#include "MMS.h"
#include "SimpleCornerBalance.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
/// SimpleCornerBalance object constructor
///
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myInput YAML input object for the simulation
SimpleCornerBalance::SimpleCornerBalance(Mesh * myMesh,\
  Materials * myMaterials,\
  YAML::Node * myInput)	      
{
  // Point to variables for mesh and input file
  mesh = myMesh;
  input = myInput;
  materials = myMaterials;

  // Size boundary value vectors to hold boundaries in each energy group
  upperBC.resize(materials->nGroups);
  lowerBC.resize(materials->nGroups);
  outerBC.resize(materials->nGroups);
  
  vector<double> inpUpperBC,inpLowerBC,inpOuterBC;

  // Check for optional inputs  
  if ((*input)["parameters"]["upperBC"]){
    inpUpperBC=(*input)["parameters"]["upperBC"].as<vector<double>>();
    
    // Check if a blanket or unique conditions are input
    if (inpUpperBC.size() == 1)
      std::fill(upperBC.begin(),upperBC.end(),inpUpperBC[0]);
    else
      upperBC = inpUpperBC;
  } 
  else{
  
    // Set default value
    std::fill(upperBC.begin(),upperBC.end(),0.0);
  } 

  if ((*input)["parameters"]["lowerBC"]){
    inpLowerBC=(*input)["parameters"]["lowerBC"].as<vector<double>>();
 
    // Check if a blanket or unique conditions are input
    if (inpLowerBC.size() == 1)
      std::fill(lowerBC.begin(),lowerBC.end(),inpLowerBC[0]);
    else
      lowerBC = inpLowerBC;

  } 
  else{
    
    // Set default value
    std::fill(lowerBC.begin(),lowerBC.end(),0.0);
  } 

  if ((*input)["parameters"]["outerBC"]){
    inpOuterBC=(*input)["parameters"]["outerBC"].as<vector<double>>(); 
    
    // Check if a blanket or unique conditions are input
    if (inpOuterBC.size() == 1)
      std::fill(outerBC.begin(),outerBC.end(),inpOuterBC[0]);
    else
      outerBC = inpOuterBC;

  } 
  else{
    
    // Set default value
    std::fill(outerBC.begin(),outerBC.end(),0.0);
  } 

};

//==============================================================================

//==============================================================================
/// Calculate the angular flux for a given source
///
/// @param [out] aFlux Angular flux solutions are stored here
/// @param [in] halfAFlux Half angle fluxes needed to solve transport equation
/// @param [in] source Source in each cell
/// @param [in] alpha Alpha in each cell
/// @param [in] energyGroup Energy group associated with this solve. Used in 
/// determining which nuclear data to use
void SimpleCornerBalance::solve(cube * aFlux,\
  cube * halfAFlux,\
  Eigen::MatrixXd * source,\
  Eigen::MatrixXd * alpha,\
  int energyGroup)
{
  for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi){

    for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu){
    
      solveAngularFlux(aFlux,halfAFlux,source,alpha,energyGroup,iXi,iMu);

    } //iMu
  } //iXi
};
//==============================================================================

//==============================================================================
void SimpleCornerBalance::solveAngularFlux(cube * aFlux,\
  cube * halfAFlux,\
  Eigen::MatrixXd * source,\
  Eigen::MatrixXd * alpha,\
  int energyGroup,int iXi,int iMu){

  double mu = mesh->quadrature[iXi].quad[iMu][1];
  
  if (mu > 0)
    solveAngularFluxPosMu(aFlux,halfAFlux,source,alpha,energyGroup,iXi,iMu);
  else
    solveAngularFluxNegMu(aFlux,halfAFlux,source,alpha,energyGroup,iXi,iMu);

}
//==============================================================================

//==============================================================================
void SimpleCornerBalance::solveAngularFluxNegMu(cube * aFlux,\
  cube * halfAFlux,\
  Eigen::MatrixXd * source,\
  Eigen::MatrixXd * alpha,\
  int energyGroup,int iXi,int iMu){


  // Index xi, mu, and weight values are stored in quadLevel object
  const int xiIndex=0,muIndex=1,etaIndex=2,weightIndex=3;

  // Temporary variables used for looping though quad set
  double xi,sigTEff,mu,alphaMinusOneHalf,alphaPlusOneHalf,weight,\
    tau,gamma,angRedistCoeff = 0,sigTEps=1E-4;
  
  xi = mesh->quadrature[iXi].quad[0][xiIndex];
  mu = mesh->quadrature[iXi].quad[iMu][muIndex];

  // Get neutron velocity in energyGroup 
  double v = materials->neutV(energyGroup);
  int zStart,rStart,zInc,borderCellZ,borderCellR,angIdx,zStartCell,\
    rStartCell;

  // Number of rows and columns for simple corner balance
  int rows = 4,cols = 4;

  // Vectors that hold the indices to be used to define upstream and
  // downstream values depending on the ordinated  
  vector<int> withinUpstreamR(2),outUpstreamR(2),withinUpstreamZ(2),\
    outUpstreamZ(2);

  // Offsets to access each corner
  Eigen::MatrixXi cornerOffset(4,2);
  Eigen::MatrixXd sigT = Eigen::MatrixXd::Zero(rows,cols);

  // Within cell leakage matrices in R and Z directions
  double kRCoeff,kZCoeff;
  Eigen::MatrixXd kR = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd kZ = Eigen::MatrixXd::Zero(rows,cols);

  // Out of cell leakage matrices in Rand Z directions
  double lRCoeff,lZCoeff;
  Eigen::MatrixXd lR = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd lZ = Eigen::MatrixXd::Zero(rows,cols);

  // Reaction matrices
  double tCoeff,RCoeff;
  Eigen::MatrixXd t = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(rows,cols);

  // A, x, and b matrices of linear system Ax=b
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(rows);
  Eigen::VectorXd x = Eigen::VectorXd::Zero(rows);
  
  // Mask matrix to select certain columns of leakage matrices.
  // Columns chosen based on upwinding scheme of simple corner
  // balance approach.
  Eigen::MatrixXd mask = Eigen::MatrixXd::Zero(rows,cols);

  // Matrix to hold downstream values
  Eigen::MatrixXd downstream = Eigen::MatrixXd::Zero(rows,cols);

  // Vector to hold upstream values
  Eigen::VectorXd upstream = Eigen::VectorXd::Zero(rows);

  // Vector to hold half-angle fluxes used in angular redistribution term 
  Eigen::VectorXd cellHalfAFlux = Eigen::VectorXd::Zero(rows);

  // Vector to hold source values in each corner 
  Eigen::VectorXd q = Eigen::VectorXd::Ones(rows);

  // Vector holding volume of each corner. Used to calculate
  // cell average value.
  Eigen::VectorXd subCellVol = Eigen::VectorXd::Zero(rows);

  // Set dirichlet boundary conditions
  double rBC,zBC;

  // MMS terms
  Eigen::VectorXd mmsQ;
  
  // Assign differencing coefficients, tau, and weight for this
  // ordinate to temporary variables 
  alphaPlusOneHalf = mesh->quadrature[iXi].alpha[iMu+1];
  alphaMinusOneHalf = mesh->quadrature[iXi].alpha[iMu];
  weight = mesh->quadrature[iXi].quad[iMu][weightIndex]; 
  tau = mesh->quadrature[iXi].tau[iMu];

  // This is the index [aFlux(:,:,angIdx)] that contains the angular 
  // flux for this ordinate
  angIdx=mesh->quadrature[iXi].ordIdx[iMu];

  // Reset corner offset
  cornerOffset.setZero();

  // Define parameters for marching from the outer radius inward
  rStart = mesh->drsCorner.size()-1;
  rStartCell = mesh->drs.size()-1;
  borderCellR = 1;
  
  // Corners whose radial boundaries are defined by values 
  // within the cell 
  withinUpstreamR = {0,3};

  // Corners whose radial boundaries are defined by values 
  // outside the cell 
  outUpstreamR = {1,2};

  // Set dirichlet boundary condition
  rBC = outerBC[energyGroup];

  // Set corner offset values defined when mu < 0
  cornerOffset(0,0) = -1;
  cornerOffset(3,0) = -1;

  // Depending on xi, define parameters for marching across the
  // axial domain	
  if (xi > 0) {

    // Marching from the bottom to the top
    zStart = 0;
    zStartCell = 0;
    zInc = 2;
    borderCellZ = -1;
    
    // Corners whose axial boundaries are defined by values
    // within the cell 
    withinUpstreamZ = {2,3};
    
    // Corners whose axial boundaries are defined by values 
    // outside the cell
    outUpstreamZ = {0,1};
    
    // Set dirichlet bc
    zBC = lowerBC[energyGroup];
  
    // Set corner offset values defined when xi > 0
    cornerOffset(2,1) = 1;
    cornerOffset(3,1) = 1;
  }
  else {			

    // Marching from the top to the bottom
    zStart = mesh->dzsCorner.size()-1;
    zStartCell = mesh->dzs.size()-1;
    zInc = -2;
    borderCellZ = 1;
  
    // Corners whose axial boundaries are defined by values
    // within the cell 
    withinUpstreamZ = {0,1};

    // Corners whose axial boundaries are defined by values 
    // outside the cell
    outUpstreamZ = {2,3};

    // Set dirichlet bc
    zBC = upperBC[energyGroup];

    // Set corner offset values defined when xi > 0
    cornerOffset(0,1) = -1;
    cornerOffset(1,1) = -1;
  }

  for (int iR = rStart, iCellR = rStartCell, countR = 0;\
    countR < mesh->drsCorner.size();\
    iR = iR - 2, --iCellR,countR = countR + 2){
    
    for (int iZ = zStart, iCellZ = zStartCell,countZ = 0;\
      countZ < mesh->dzsCorner.size();\
      iZ = iZ + zInc,iCellZ = iCellZ + zInc/2,countZ = countZ + 2){
     
      // Set source in each corner
      for (int iCorner = 0; iCorner < 4; ++iCorner){
        q(iCorner) = (*source)(iZ+cornerOffset(iCorner,1),\
        iR+cornerOffset(iCorner,0));
      }

      for (int iSig = 0; iSig < sigT.cols(); ++iSig){
        sigTEff = materials->sigT(iZ+cornerOffset(iSig,1),\
          iR+cornerOffset(iSig,0),energyGroup)\
          +(*alpha)(iZ+cornerOffset(iSig,1),iR+cornerOffset(iSig,0))/v; 

        if (sigTEff > sigTEps)
          sigT(iSig,iSig) = sigTEff;
        else
          sigT(iSig,iSig) = sigTEps;
      }


      // Calculate the the ratio of the inner to outer radius
      // for this cell
      gamma = mesh->rEdge(iCellR)/mesh->rEdge(iCellR+1);
      
      // Calculate radial within cell leakage matrix
      kRCoeff = mesh->dzs(iCellZ)*mesh->rEdge(iCellR+1)/8.0;
      kR=calckR(gamma);
      kR=kRCoeff*kR;
      
      // Calculate axial within cell leakage matrix
      kZCoeff = mesh->drs(iCellR)*mesh->rEdge(iCellR+1)/16.0;
      kZ=calckZ(gamma);
      kZ=kZCoeff*kZ;
      
      // Calculate radial out of cell leakage matrix
      lRCoeff = mesh->dzs(iCellZ)*mesh->rEdge(iCellR+1)/2.0;
      lR=calclR(gamma);
      lR=lRCoeff*lR;
                              
      // Calculate axial out of cell leakage matrix
      lZCoeff = mesh->drs(iCellR)*mesh->rEdge(iCellR+1)/8.0;
      lZ=calclZ(gamma);
      lZ=lZCoeff*lZ;
      
      // Calculate first collision matrix
      tCoeff = mesh->drs(iCellR)*mesh->dzs(iCellZ)\
        *mesh->rEdge(iCellR+1)/16.0;
      t=calct(gamma);
      t=tCoeff*t;
                              
      // Calculate angular redistribution matrix
      RCoeff = mesh->drs(iCellR)*mesh->dzs(iCellZ)/4.0;
      R=calcR(gamma);
      R=RCoeff*R;
      angRedistCoeff = alphaPlusOneHalf/(weight*tau); 
      

      cout << "initial A: " << endl;   
      cout << A << endl; 

      // Calculate A considering within cell leakage, collision,
      // and angular redistribution
      A = mu*kR+xi*kZ+sigT*t+angRedistCoeff*R;
  
      cout << "A: " << endl;   
      cout << A << endl; 
      
      // Consider radial boundary values defined in this cell
      mask.setIdentity();
      for (int iCol = 0; iCol < outUpstreamR.size(); ++iCol){
        mask(outUpstreamR[iCol],outUpstreamR[iCol])=0;
      }
      downstream = mu*lR*mask;	
      
      // Consider axial boundary values defined in this cell
      mask.setIdentity();
      for (int iCol = 0; iCol < outUpstreamZ.size(); ++iCol){
        mask(outUpstreamZ[iCol],outUpstreamZ[iCol])=0;
      }
      downstream = downstream + xi*lZ*mask;	
      
      A = A + downstream;
      
      cout << "boundary values A: " << endl;   
      cout << A << endl; 
      
      cout << "initial b: " << endl;   
      cout << b << endl; 
      
      // Form b matrix
      b = t*q;
      
      cout << "source b: " << endl;   
      cout << b << endl; 


      // Consider contribution of angular redistribution term
      // calculated with known values
      angRedistCoeff = ((alphaPlusOneHalf/tau)*(tau - 1.0)\
        - alphaMinusOneHalf)/weight;
      cellHalfAFlux.setOnes();
      cellHalfAFlux=(*halfAFlux)(iZ,iR,iXi)*cellHalfAFlux;
      
      // Read corner half angle fluxes into a vector
      for (int iCorner = 0; iCorner < 4; ++iCorner){
        cellHalfAFlux(iCorner) = (*halfAFlux)(iZ+cornerOffset(iCorner,1),\
        iR+cornerOffset(iCorner,0),iXi);
      }

      b = b - angRedistCoeff*R*cellHalfAFlux;
      
      cout << "ang redist b: " << endl;   
      cout << b << endl; 
      
      // Consider radial boundary values defined in other cells 
      // or by BCs
      if (iR!=rStart){

        upstream =\
        mu*(*aFlux)(iZ+cornerOffset(outUpstreamR[0],1),\
        iR+borderCellR,angIdx)*lR.col(outUpstreamR[0])+\
        mu*(*aFlux)(iZ+cornerOffset(outUpstreamR[1],1),\
        iR+borderCellR,angIdx)*lR.col(outUpstreamR[1]); 

        b = b - upstream;

      } else {
        upstream = mu*rBC\
        *(lR.col(outUpstreamR[0])+lR.col(outUpstreamR[1]));
        b = b - upstream;
      }
      
      cout << "radial boundary b: " << endl;   
      cout << b << endl; 
      
      // Consider axial boundary values defined in other cells 
      // or by BCs
      if (iZ!=zStart){

        upstream =\
        xi*(*aFlux)(iZ+borderCellZ,iR+cornerOffset(outUpstreamZ[0],0),\
        angIdx)*lZ.col(outUpstreamZ[0])+\
        xi*(*aFlux)(iZ+borderCellZ,iR+cornerOffset(outUpstreamZ[1],0),\
        angIdx)*lZ.col(outUpstreamZ[1]);

        b = b - upstream;

      } else{
        upstream = xi*zBC\
        *(lZ.col(outUpstreamZ[0])+lZ.col(outUpstreamZ[1]));
        b = b - upstream;
      }
      
      cout << "axial boundary b: " << endl;   
      cout << b << endl; 
      
      subCellVol = calcSubCellVol(iCellZ,iCellR);	

      // Solve for angular fluxes in each corner
      x = A.partialPivLu().solve(b);

      cout <<"x: "<< endl;
      cout << x << endl;

      // Take average of corner values to get angular flux 
      // for this cell
      subCellVol = calcSubCellVol(iCellZ,iCellR);	
      
      for (int iCorner = 0; iCorner < 4; ++iCorner){
        (*aFlux)(iZ+cornerOffset(iCorner,1),\
        iR+cornerOffset(iCorner,0),angIdx) = x(iCorner); 
      }      
      // Use weighted diamond difference to calculate next half
      // angle flux used for next value of mu in this quadrature 
      // level
      for (int iCorner = 0; iCorner < 4; ++iCorner){
        (*halfAFlux)(iZ+cornerOffset(iCorner,1),\
        iR+cornerOffset(iCorner,0),iXi) = ((*aFlux)\
        (iZ+cornerOffset(iCorner,1),iR+cornerOffset(iCorner,0),angIdx)\
        +(tau-1.0)*(*halfAFlux)(iZ+cornerOffset(iCorner,1),\
        iR+cornerOffset(iCorner,0),iXi))/tau; 
      }      
    } //iR
  } //iZ
}
//==============================================================================

//==============================================================================
void SimpleCornerBalance::solveAngularFluxPosMu(cube * aFlux,\
  cube * halfAFlux,\
  Eigen::MatrixXd * source,\
  Eigen::MatrixXd * alpha,\
  int energyGroup,int iXi,int iMu){


  // Index xi, mu, and weight values are stored in quadLevel object
  const int xiIndex=0,muIndex=1,etaIndex=2,weightIndex=3;

  // Temporary variables used for looping though quad set
  double xi,sigTEff,mu,alphaMinusOneHalf,alphaPlusOneHalf,weight,\
    tau,gamma,angRedistCoeff = 0,sigTEps=1E-4;
  
  xi = mesh->quadrature[iXi].quad[0][xiIndex];
  mu = mesh->quadrature[iXi].quad[iMu][muIndex];

  // Get neutron velocity in energyGroup 
  double v = materials->neutV(energyGroup);
  int numPs,numQs,reflectedP,reflectedQ,reflectedAngIdx,zStart,rStart,\
    zInc,borderCellZ,borderCellR,angIdx,zStartCell,rStartCell;

  // Number of rows and columns for simple corner balance
  int rows = 4,cols = 4;

  // Vectors that hold the indices to be used to define upstream and
  // downstream values depending on the ordinated  
  vector<int> withinUpstreamR(2),outUpstreamR(2),withinUpstreamZ(2),\
    outUpstreamZ(2);

  // Offsets to access each corner
  Eigen::MatrixXi cornerOffset(4,2);
  Eigen::MatrixXd sigT = Eigen::MatrixXd::Zero(rows,cols);

  // Within cell leakage matrices in R and Z directions
  double kRCoeff,kZCoeff;
  Eigen::MatrixXd kR = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd kZ = Eigen::MatrixXd::Zero(rows,cols);

  // Out of cell leakage matrices in Rand Z directions
  double lRCoeff,lZCoeff;
  Eigen::MatrixXd lR = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd lZ = Eigen::MatrixXd::Zero(rows,cols);

  // Reaction matrices
  double tCoeff,RCoeff;
  Eigen::MatrixXd t = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(rows,cols);

  // A, x, and b matrices of linear system Ax=b
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(rows);
  Eigen::VectorXd x = Eigen::VectorXd::Zero(rows);
  
  // Mask matrix to select certain columns of leakage matrices.
  // Columns chosen based on upwinding scheme of simple corner
  // balance approach.
  Eigen::MatrixXd mask = Eigen::MatrixXd::Zero(rows,cols);

  // Matrix to hold downstream values
  Eigen::MatrixXd downstream = Eigen::MatrixXd::Zero(rows,cols);

  // Vector to hold upstream values
  Eigen::VectorXd upstream = Eigen::VectorXd::Zero(rows);

  // Vector to hold half-angle fluxes used in angular redistribution term 
  Eigen::VectorXd cellHalfAFlux = Eigen::VectorXd::Zero(rows);

  // Vector to hold source values in each corner 
  Eigen::VectorXd q = Eigen::VectorXd::Ones(rows);

  // Vector holding volume of each corner. Used to calculate
  // cell average value.
  Eigen::VectorXd subCellVol = Eigen::VectorXd::Zero(rows);

  // Set dirichlet boundary conditions
  double rBC,zBC;

  // MMS terms
  Eigen::VectorXd mmsQ;
      
  // Assign differencing coefficients, tau, and weight for this
  // ordinate to temporary variables
  alphaPlusOneHalf = mesh->quadrature[iXi].alpha[iMu+1];
  alphaMinusOneHalf = mesh->quadrature[iXi].alpha[iMu];
  weight = mesh->quadrature[iXi].quad[iMu][weightIndex]; 
  tau = mesh->quadrature[iXi].tau[iMu];

  // This is the index [aFlux(:,:,angIdx)] that contains the angular
  // flux for this ordinate
  angIdx= mesh->quadrature[iXi].ordIdx[iMu];

  // Reset corner offset
  cornerOffset.setZero();
  
  // Define parameters for marching from the outer radius inward
  rStart = 0;
  rStartCell = 0;
  borderCellR = -1;

  // Corners whose radial boundaries are defined by values
  // within the cell
  withinUpstreamR = {1,2};

  // Corners whose radial boundaries are defined by values
  // outside the cell
  outUpstreamR = {0,3};

  // Set corner offset values defined when mu < 0
  cornerOffset(1,0) = 1;
  cornerOffset(2,0) = 1;

  // Depending on xi, define parameters for marching across the
  // axial domain
  if (xi > 0) {

    // Marching from the bottom to the top
    zStart = 0;
    zStartCell = 0;
    zInc = 2;
    borderCellZ = -1;

    // Corners whose axial boundaries are defined by values
    // within the cell
    withinUpstreamZ = {2,3};

    // Corners whose axial boundaries are defined by values
    // outside the cell
    outUpstreamZ = {0,1};
  
    // Set dirichlet boundary condition
    zBC = lowerBC[energyGroup];

    // Set corner offset values defined when xi > 0
    cornerOffset(2,1) = 1;
    cornerOffset(3,1) = 1;
  }
  else {			

    // Marching from the top to the bottom
    zStart = mesh->dzsCorner.size()-1;
    zStartCell = mesh->dzs.size()-1;
    zInc = -2;
    borderCellZ = 1;
    
    // Corners whose axial boundaries are defined by values
    // within the cell
    withinUpstreamZ = {0,1};

    // Corners whose axial boundaries are defined by values
    // outside the cell
    outUpstreamZ = {2,3};
    
    // Set dirichlet boundary condition
    zBC = upperBC[energyGroup];

    // Set corner offset values defined when xi < 0
    cornerOffset(0,1) = -1;
    cornerOffset(1,1) = -1;
  }

  for (int iR = rStart, iCellR = rStartCell, countR = 0;\
    countR < mesh->drsCorner.size();\
    iR = iR + 2, ++iCellR, countR = countR + 2){
   
     
    for (int iZ = zStart, iCellZ = zStartCell,countZ = 0;\
      countZ < mesh->dzsCorner.size();\
      iZ = iZ + zInc, iCellZ = iCellZ + zInc/2,countZ = countZ + 2){

      // Set source in each corner
      for (int iCorner = 0; iCorner < 4; ++iCorner){
        q(iCorner) = (*source)(iZ+cornerOffset(iCorner,1),\
        iR+cornerOffset(iCorner,0));
      }

      for (int iSig = 0; iSig < sigT.cols(); ++iSig){
        sigTEff = materials->sigT(iZ+cornerOffset(iSig,1),\
          iR+cornerOffset(iSig,0),energyGroup)\
          +(*alpha)(iZ+cornerOffset(iSig,1),iR+cornerOffset(iSig,0))/v; 

        if (sigTEff > sigTEps)
          sigT(iSig,iSig) = sigTEff;
        else
          sigT(iSig,iSig) = sigTEps;
      }

      // Calculate the the ratio of the inner to outer radius
      // for this cell
      gamma = mesh->rEdge(iCellR)/mesh->rEdge(iCellR+1);
      
      // Calculate radial within cell leakage matrix
      kRCoeff = mesh->dzs(iCellZ)*mesh->rEdge(iCellR+1)/8.0;
      kR=calckR(gamma);
      kR=kRCoeff*kR;
      
      // Calculate axial within cell leakage matrix
      kZCoeff = mesh->drs(iCellR)*mesh->rEdge(iCellR+1)/16.0;
      kZ=calckZ(gamma);
      kZ=kZCoeff*kZ;
      
      // Calculate radial out of cell leakage matrix
      lRCoeff = mesh->dzs(iCellZ)*mesh->rEdge(iCellR+1)/2.0;
      lR=calclR(gamma);
      lR=lRCoeff*lR;
                              
      // Calculate axial out of cell leakage matrix
      lZCoeff = mesh->drs(iCellR)*mesh->rEdge(iCellR+1)/8.0;
      lZ=calclZ(gamma);
      lZ=lZCoeff*lZ;
      
      // Calculate first collision matrix
      tCoeff = mesh->drs(iCellR)*mesh->dzs(iCellZ)\
        *mesh->rEdge(iCellR+1)/16.0;
      t=calct(gamma);
      t=tCoeff*t;
                              
      // Calculate second collision matrix
      RCoeff = mesh->drs(iCellR)*mesh->dzs(iCellZ)/4.0;
      R=calcR(gamma);
      R=RCoeff*R;
      angRedistCoeff = alphaPlusOneHalf/(weight*tau); 

      // Calculate A considering within cell leakage, collision,
      // and angular redistribution
      A = mu*kR+xi*kZ+sigT*t+angRedistCoeff*R;

      // Consider radial boundary values defined in this cell
      mask.setIdentity();
      for (int iCol = 0; iCol < outUpstreamR.size(); ++iCol){
        mask(outUpstreamR[iCol],outUpstreamR[iCol])=0;
      }
      downstream = mu*lR*mask;	
      
      // Consider axial boundary values defined in this cell
      mask.setIdentity();
      for (int iCol = 0; iCol < outUpstreamZ.size(); ++iCol){
        mask(outUpstreamZ[iCol],outUpstreamZ[iCol])=0;
      }
      downstream = downstream + xi*lZ*mask;	
      
      A = A + downstream;

      // Form b matrix
      b = t*q;
    
      // Consider contribution of angular redistribution term 
      // calculated with known values
      angRedistCoeff = ((alphaPlusOneHalf/tau)*(tau - 1.0)\
        - alphaMinusOneHalf)/weight;
      cellHalfAFlux.setOnes();
      cellHalfAFlux=(*halfAFlux)(iZ,iR,iXi)*cellHalfAFlux;

      // Read corner half angle fluxes into a vector
      for (int iCorner = 0; iCorner < 4; ++iCorner){
        cellHalfAFlux(iCorner) = (*halfAFlux)(iZ+cornerOffset(iCorner,1),\
        iR+cornerOffset(iCorner,0),iXi);
      }

      b = b - angRedistCoeff*R*cellHalfAFlux;

      // Consider radial boundary values in other cells or BCs
      if (iR==rStart){
        // Apply reflecting boundary condition

        // Calculate indices of reflected ordinate
        numPs = mesh->quadrature.size();
        numQs = mesh->quadrature[iXi].nOrd;
        reflectedP= numPs-iXi-1;
        reflectedQ= numQs-iMu-1;

        // Get angular index of reflect flux
        reflectedAngIdx = mesh->quadrature[iXi].ordIdx[reflectedQ]; 
        
        upstream =\
        mu*(*aFlux)(iZ+cornerOffset(outUpstreamR[0],1),\
        iR,reflectedAngIdx)*lR.col(outUpstreamR[0])+
        mu*(*aFlux)(iZ+cornerOffset(outUpstreamR[1],1),\
        iR,reflectedAngIdx)*lR.col(outUpstreamR[1]);

        b = b - upstream;
      }

      // Consider radial boundary values defined in other cells
      // or by BCs
      else if (iR!=rStart){
        upstream =\
        mu*(*aFlux)(iZ+cornerOffset(outUpstreamR[0],1),\
        iR+borderCellR,angIdx)*lR.col(outUpstreamR[0])+
        mu*(*aFlux)(iZ+cornerOffset(outUpstreamR[1],1),\
        iR+borderCellR,angIdx)*lR.col(outUpstreamR[1]);
        b = b - upstream;
      
      }

      // Consider axial boundary values defined in other cells
      // or by BCs
      if (iZ!=zStart){
        upstream =\
        xi*(*aFlux)(iZ+borderCellZ,iR+cornerOffset(outUpstreamZ[0],0),\
        angIdx)*lZ.col(outUpstreamZ[0])+\
        xi*(*aFlux)(iZ+borderCellZ,iR+cornerOffset(outUpstreamZ[1],0),\
        angIdx)*lZ.col(outUpstreamZ[1]);

        b = b - upstream;
      } else{
        upstream = xi*zBC\
        *(lZ.col(outUpstreamZ[0])+lZ.col(outUpstreamZ[1]));
        b = b - upstream;
      }
                  
      subCellVol = calcSubCellVol(iCellZ,iCellR);	

      // Solve for angular fluxes in each corner 
      x = A.partialPivLu().solve(b);

      // Take average of corner values to get angular flux
      // for this cell
      subCellVol = calcSubCellVol(iCellZ,iCellR);

      for (int iCorner = 0; iCorner < 4; ++iCorner){
        (*aFlux)(iZ+cornerOffset(iCorner,1),\
        iR+cornerOffset(iCorner,0),angIdx) = x(iCorner); 
      }      
    
      // Use weighted diamond difference to calculate next half
      // angle flux used for next nvalue of mu in this quadrature 
      // level
      for (int iCorner = 0; iCorner < 4; ++iCorner){
        (*halfAFlux)(iZ+cornerOffset(iCorner,1),\
        iR+cornerOffset(iCorner,0),iXi) = ((*aFlux)\
        (iZ+cornerOffset(iCorner,1),iR+cornerOffset(iCorner,0),angIdx)\
        +(tau-1.0)*(*halfAFlux)(iZ+cornerOffset(iCorner,1),\
        iR+cornerOffset(iCorner,0),iXi))/tau; 
      }      
    } //iR
  } //iZ
}
//==============================================================================




//==============================================================================
/// Calculate within cell radial leakage matrix
///
/// @param [in] myGamma Ratio of inner radius to outer radius in a cell 
/// @param [out] kR Within cell radial leakage matrix
Eigen::MatrixXd SimpleCornerBalance::calckR(double myGamma){
  double a = (1+myGamma);
  double b = -(1+myGamma);
  Eigen::MatrixXd kR = Eigen::MatrixXd::Zero(4,4);

  kR(0,0) = a; kR(0,1) = a;
  kR(1,0) = b; kR(1,1) = b;
  kR(2,2) = b; kR(2,3) = b;
  kR(3,2) = a; kR(3,3) = a;
           
  return kR;
}

//==============================================================================


//==============================================================================
/// Calculate within cell axial leakage matrix
///
/// @param [in] myGamma Ratio of inner radius to outer radius in a cell 
/// @param [out] kZ Within cell axial leakage matrix
Eigen::MatrixXd SimpleCornerBalance::calckZ(double myGamma){
  double a = 1+3*myGamma;
  double b = 3+myGamma;
  Eigen::MatrixXd kZ = Eigen::MatrixXd::Zero(4,4);

  kZ(0,0) = a; kZ(0,3) = a;
  kZ(1,1) = b; kZ(1,2) = b;
  kZ(2,1) = -b; kZ(2,2) = -b;
  kZ(3,0) = -a; kZ(3,3) = -a;

  return kZ;
}

//==============================================================================
/// Calculate out of cell radial leakage matrix
///
/// @param [in] myGamma Ratio of inner radius to outer radius in a cell 
/// @param [out] lR Out of cell radial leakage matrix
Eigen::MatrixXd SimpleCornerBalance::calclR(double myGamma){
  double a = -myGamma;
  double b = 1;
  Eigen::MatrixXd lR = Eigen::MatrixXd::Zero(4,4);

  lR(0,0) = a; 
  lR(1,1) = b; 
  lR(2,2) = b; 
  lR(3,3) = a; 
           
  return lR;
}
//==============================================================================

//==============================================================================
/// Calculate out of cell axial leakage matrix
///
/// @param [in] myGamma Ratio of inner radius to outer radius in a cell 
/// @param [out] lZ Out of cell axial leakage matrix
Eigen::MatrixXd SimpleCornerBalance::calclZ(double myGamma){
  double a = 1+3*myGamma;
  double b = 3+myGamma;
  Eigen::MatrixXd lZ = Eigen::MatrixXd::Zero(4,4);

  lZ(0,0) = -a; 
  lZ(1,1) = -b; 
  lZ(2,2) = b; 
  lZ(3,3) = a; 
           
  return lZ;
}
//==============================================================================

//==============================================================================
/// Calculate collision matrix
///
/// @param [in] myGamma Ratio of inner radius to outer radius in a cell 
/// @param [out] t Collision matrix
Eigen::MatrixXd SimpleCornerBalance::calct(double myGamma){
  double a = 1+3*myGamma;
  double b = 3+myGamma;
  Eigen::MatrixXd t = Eigen::MatrixXd::Zero(4,4);

  t(0,0) = a; 
  t(1,1) = b; 
  t(2,2) = b; 
  t(3,3) = a; 
           
  return t;
}
//==============================================================================

//==============================================================================
/// Calculate angular redistribution matrix
///
/// @param [in] myGamma Ratio of inner radius to outer radius in a cell 
/// @param [out] R Angular redistribution matrix
Eigen::MatrixXd SimpleCornerBalance::calcR(double myGamma){
  double a = 1;
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(4,4);

  R(0,0) = a; 
  R(1,1) = a; 
  R(2,2) = a;
  R(3,3) = a; 
           
  return R;
}
//==============================================================================

//==============================================================================
/// Calculate volumes of subcell regions
///
/// @param [out] subCellVol Contains the volume in each cell
Eigen::VectorXd SimpleCornerBalance::calcSubCellVol(int myiZ, int myiR){
  Eigen::VectorXd subCellVol = Eigen::VectorXd::Zero(4);
          
  subCellVol(0) = (mesh->dzs(myiZ)/2)*(pow(mesh->rCent(myiR),2)-\
    pow(mesh->rEdge(myiR),2))/2;

  subCellVol(1) = (mesh->dzs(myiZ)/2)*(pow(mesh->rEdge(myiR+1),2)-\
    pow(mesh->rCent(myiR),2))/2;

  subCellVol(2) = (mesh->dzs(myiZ)/2)*(pow(mesh->rEdge(myiR+1),2)-\
    pow(mesh->rCent(myiR),2))/2;

  subCellVol(3) = (mesh->dzs(myiZ)/2)*(pow(mesh->rCent(myiR),2)-\
    pow(mesh->rEdge(myiR),2))/2;

  return subCellVol;
}
//==============================================================================


