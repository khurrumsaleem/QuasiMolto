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
#include "SimpleCornerBalance.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
//! SimpleCornerBalance object constructor

SimpleCornerBalance::SimpleCornerBalance(Mesh * myMesh,\
  Materials * myMaterials,\
  YAML::Node * myInput)	      
{
// Point to variables for mesh and input file
  mesh = myMesh;
  input = myInput;
  materials = myMaterials;
};

//==============================================================================

//==============================================================================
//! solve Calculate the angular flux for a given source

// Solves using simple corner balance on the neutron transport equation 
// in RZ geometry. 

void SimpleCornerBalance::solve(cube * aFlux,\
  cube * halfAFlux,\
  Eigen::MatrixXd * source,\
  int energyGroup)
{
  // index xi, mu, and weight values are stored in quadLevel object
  const int xiIndex=0,muIndex=1,weightIndex=3;

  // temporary variables used for looping though quad set
  double xi,sigT,mu,alphaMinusOneHalf,alphaPlusOneHalf,weight,\
    tau,gamma,angRedistCoeff = 0;
  int numPs,numQs,reflectedP,reflectedQ,reflectedAngIdx,zStart,rStart,\
    zEnd,zInc,borderCellZ,borderCellR,angIdx;

  // number of rows and columns for simple corner balance
  int rows = 4,cols = 4;

  // vectors that hold the indices to be used to define upstream and
  // downstream values depending on the ordinated  
  vector<int> withinUpstreamR(2),outUpstreamR(2),withinUpstreamZ(2),\
    outUpstreamZ(2);

  // within cell leakage matrices in R and Z directions
  double kRCoeff,kZCoeff;
  Eigen::MatrixXd kR = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd kZ = Eigen::MatrixXd::Zero(rows,cols);

  // out of cell leakage matrices in Rand Z directions
  double lRCoeff,lZCoeff;
  Eigen::MatrixXd lR = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd lZ = Eigen::MatrixXd::Zero(rows,cols);

  // reaction matrices
  double tCoeff,RCoeff;
  Eigen::MatrixXd t = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(rows,cols);

  // A, x, and b matrices of linear system Ax=b
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(rows);
  Eigen::VectorXd x = Eigen::VectorXd::Zero(rows);
  
  // mask matrix to select certain columns of leakage matrices.
  // columns chosen based on upwinding scheme of simple corner
  // balance approach.
  Eigen::MatrixXd mask = Eigen::MatrixXd::Zero(rows,cols);

  // matrix to hold downstream values
  Eigen::MatrixXd downstream = Eigen::MatrixXd::Zero(rows,cols);

  // vector to hold upstream values
  Eigen::VectorXd upstream = Eigen::VectorXd::Zero(rows);

  // vector to hold half-angle fluxes used in angular redistribution term 
  Eigen::VectorXd cellHalfAFlux = Eigen::VectorXd::Zero(rows);

  // vector to hold source values in each corner 
  Eigen::VectorXd q = Eigen::VectorXd::Ones(rows);

  // vector holding volume of each corner. Used to calculate
  // cell average value.
  Eigen::VectorXd subCellVol = Eigen::VectorXd::Zero(rows);

  for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi){

    // get xi for this quadrature level
    xi = mesh->quadrature[iXi].quad[0][xiIndex];

    for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu){

      // get xi for this ordinate 
      mu = mesh->quadrature[iXi].quad[iMu][muIndex];

      // assign differencing coefficients, tau, and weight for this
      // ordinate to temporary variables 
      alphaPlusOneHalf = mesh->quadrature[iXi].alpha[iMu+1];
      alphaMinusOneHalf = mesh->quadrature[iXi].alpha[iMu];
      weight = mesh->quadrature[iXi].quad[iMu][weightIndex]; 
      tau = mesh->quadrature[iXi].tau[iMu];
      
      // this is the index [aFlux(:,:,angIdx)] that contains the angular 
      // flux for this ordinate
      angIdx= mesh->quadrature[iXi].ordIdx[iMu];

      // do cases where mu < 0 first, as those solutions are needed to 
      // to define the reflecting boundary condition at Z=0; 
      if (mu < 0 ){
        
        // define parameters for marching from the outer radius inward
        rStart = mesh->drs.size()-1;
        borderCellR = 1;
        
        // corners whose radial boundaries are defined by values 
        // within the cell 
        withinUpstreamR = {0,3};

        // corners whose radial boundaries are defined by values 
        // outside the cell 
        outUpstreamR = {1,2};

        // depending on xi, define parameters for marching across the
        // axial domain	
        if (xi > 0) {
 
          // marching from the top to the bottom
          zStart = mesh->dzs.size()-1;
          zEnd = 0;
          zInc = -1;
          borderCellZ = 1;
          withinUpstreamZ = {2,3};
          outUpstreamZ = {0,1};
        }
        else {			

          // marching from the bottom to the top
          zStart = 0;
          zEnd = mesh->dzs.size();
          zInc = 1;
          borderCellZ = -1;
        
          // corners whose axial boundaries are defined by values
          // within the cell 
          withinUpstreamZ = {0,1};

          // corners whose axial boundaries are defined by values 
          // outside the cell
          outUpstreamZ = {2,3};
        }

        for (int iR = rStart, countR = 0; countR < mesh->drs.size();\
          --iR, ++countR){

          for (int iZ = zStart, countZ = 0; \
            countZ < mesh->dzs.size(); iZ = iZ + zInc, ++countZ){
            
            // calculate source in each corner
            q.setOnes();
            q = (*source)(iZ,iR)*q;

            // get the total cross section in this cell
            sigT = materials->sigT(iZ,iR,energyGroup);

            // calculate the the ratio of the inner to outer radius
            // for this cell
            gamma = mesh->rEdge(iR)/mesh->rEdge(iR+1);
            
            // calculate radial within cell leakage matrix
            kRCoeff = mesh->dzs(iZ)*mesh->rEdge(iR+1)/8.0;
            kR=calckR(gamma);
            kR=kRCoeff*kR;
            
            // calculate axial within cell leakage matrix
            kZCoeff = mesh->drs(iR)*mesh->rEdge(iR+1)/16.0;
            kZ=calckZ(gamma);
            kZ=kZCoeff*kZ;
            
            // calculate radial out of cell leakage matrix
            lRCoeff = mesh->dzs(iZ)*mesh->rEdge(iR+1)/2.0;
            lR=calclR(gamma);
            lR=lRCoeff*lR;
                                    
            // calculate axial out of cell leakage matrix
            lZCoeff = mesh->drs(iR)*mesh->rEdge(iR+1)/8.0;
            lZ=calclZ(gamma);
            lZ=lZCoeff*lZ;
            
            // calculate first collision matrix
            tCoeff = mesh->drs(iR)*mesh->dzs(iZ)*mesh->rEdge(iR+1)/16.0;
            t=calct(gamma);
            t=tCoeff*t;
                                    
            // calculate angular redistribution matrix
            RCoeff = mesh->drs(iR)*mesh->dzs(iZ)/4.0;
            R=calcR(gamma);
            R=RCoeff*R;
            angRedistCoeff = alphaPlusOneHalf/(weight*tau); 

            // calculate A considering within cell leakage, collision,
            // and angular redistribution
            A = mu*kR+xi*kZ+sigT*t+angRedistCoeff*R;

            // consider radial boundary values defined in this cell
            mask.setIdentity();
            for (int iCol = 0; iCol < outUpstreamR.size(); ++iCol){
              mask(outUpstreamR[iCol],outUpstreamR[iCol])=0;
            }
            downstream = mu*lR*mask;	
            
            // consider axial boundary values defined in this cell
            mask.setIdentity();
            for (int iCol = 0; iCol < outUpstreamZ.size(); ++iCol){
              mask(outUpstreamZ[iCol],outUpstreamZ[iCol])=0;
            }
            downstream = downstream + xi*lZ*mask;	
            
            A = A + downstream;

            // form b matrix
            b = t*q;

            // consider contribution of angular redistribution term
            // calculated with known values
            angRedistCoeff = ((alphaPlusOneHalf/tau)*(tau - 1.0)\
              - alphaMinusOneHalf)/weight;
            cellHalfAFlux.setOnes();
            cellHalfAFlux=(*halfAFlux)(iZ,iR,iXi)*cellHalfAFlux;
            b = b - angRedistCoeff*R*cellHalfAFlux;

            // consider radial boundary values defined in other cells 
            // or by BCs
            if (iR!=rStart){
              upstream = mu*(*aFlux)(iZ,iR+borderCellR,angIdx)\
              *(lR.col(outUpstreamR[0])+lR.col(outUpstreamR[1]));
              b = b - upstream;
            }

            // consider axial boundary values defined in other cells 
            // or by BCs
            if (iZ!=zStart){
              upstream = xi*(*aFlux)(iZ+borderCellZ,iR,angIdx)\
              *(lZ.col(outUpstreamZ[0])+lZ.col(outUpstreamZ[1]));
              b = b - upstream;
            }

            // solve for angular fluxes in each corner
            x = A.lu().solve(b);

            // take average of corner values to get angular flux 
            // for this cell
            subCellVol = calcSubCellVol(iZ,iR);	
            (*aFlux)(iZ,iR,angIdx) = x.dot(subCellVol)/subCellVol.sum();

            // use weighted diamond difference to calculate next half
            // angle flux used for next value of mu in this quadrature 
            // level
            (*halfAFlux)(iZ,iR,iXi) = ((*aFlux)(iZ,iR,angIdx)\
            +(tau-1.0)*(*halfAFlux)(iZ,iR,iXi))/tau; 

          } //iR
        } //iZ
      } // if mu < 0
    } //iMu
  } //iXi
  
  // Repeat for mu > 0
  for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi){

    // get xi for this quadrature level
    xi = mesh->quadrature[iXi].quad[0][xiIndex];

    for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu){
    
      // get xi for this ordinate 
      mu = mesh->quadrature[iXi].quad[iMu][muIndex]; 

      // assign differencing coefficients, tau, and weight for this
      // ordinate to temporary variables
      alphaPlusOneHalf = mesh->quadrature[iXi].alpha[iMu+1];
      alphaMinusOneHalf = mesh->quadrature[iXi].alpha[iMu];
      weight = mesh->quadrature[iXi].quad[iMu][weightIndex]; 
      tau = mesh->quadrature[iXi].tau[iMu];

      // this is the index [aFlux(:,:,angIdx)] that contains the angular
      // flux for this ordinate
      angIdx= mesh->quadrature[iXi].ordIdx[iMu];

      // do cases where mu > 0, as the boundary values are defined now
      if (mu > 0 ){

        // define parameters for marching from the outer radius inward
        rStart = 0;
        borderCellR = -1;

        // corners whose radial boundaries are defined by values
        // within the cell
        withinUpstreamR = {1,2};

        // corners whose radial boundaries are defined by values
        // outside the cell
        outUpstreamR = {0,3};

        // depending on xi, define parameters for marching across the
        // axial domain
        if (xi > 0) {

          // marhcing from the top to the bottom
          zStart = mesh->dzs.size()-1;
          zEnd = 0;
          zInc = -1;
          borderCellZ = 1;

          // corners whose axial boundaries are defined by values
          // within the cell
          withinUpstreamZ = {2,3};

          // corners whose axial boundaries are defined by values
          // outside the cell
          outUpstreamZ = {0,1};
        }
        else {			

          // marching from the bottom to the top
          zStart = 0;
          zEnd = mesh->dzs.size();
          zInc = 1;
          borderCellZ = -1;
          
          // corners whose axial boundaries are defined by values
          // within the cell
          withinUpstreamZ = {0,1};

          // corners whose axial boundaries are defined by values
          // outside the cell
          outUpstreamZ = {2,3};
        }
        for (int iR = rStart, countR = 0; countR < mesh->drs.size();\
           ++iR, ++countR){
          
          for (int iZ = zStart, countZ = 0; \
            countZ < mesh->dzs.size(); iZ = iZ + zInc, ++countZ){

            // calculate source in each corner
            q.setOnes();
            q = (*source)(iZ,iR)*q;

            // get the total cross section in this cell
            sigT = materials->sigT(iZ,iR,energyGroup);

            // calculate the the ratio of the inner to outer radius
            // for this cell
            gamma = mesh->rEdge(iR)/mesh->rEdge(iR+1);
            
            // calculate radial within cell leakage matrix
            kRCoeff = mesh->dzs(iZ)*mesh->rEdge(iR+1)/8.0;
            kR=calckR(gamma);
            kR=kRCoeff*kR;
            
            // calculate axial within cell leakage matrix
            kZCoeff = mesh->drs(iR)*mesh->rEdge(iR+1)/16.0;
            kZ=calckZ(gamma);
            kZ=kZCoeff*kZ;
            
            // calculate radial out of cell leakage matrix
            lRCoeff = mesh->dzs(iZ)*mesh->rEdge(iR+1)/2.0;
            lR=calclR(gamma);
            lR=lRCoeff*lR;
                                    
            // calculate axial out of cell leakage matrix
            lZCoeff = mesh->drs(iR)*mesh->rEdge(iR+1)/8.0;
            lZ=calclZ(gamma);
            lZ=lZCoeff*lZ;
            
            // calculate first collision matrix
            tCoeff = mesh->drs(iR)*mesh->dzs(iZ)*mesh->rEdge(iR+1)/16.0;
            t=calct(gamma);
            t=tCoeff*t;
                                    
            // calculate second collision matrix
            RCoeff = mesh->drs(iR)*mesh->dzs(iZ)/4.0;
            R=calcR(gamma);
            R=RCoeff*R;
            angRedistCoeff = alphaPlusOneHalf/(weight*tau); 

            // calculate A considering within cell leakage, collision,
            // and angular redistribution
            A = mu*kR+xi*kZ+sigT*t+angRedistCoeff*R;

            // consider radial boundary values defined in this cell
            mask.setIdentity();
            for (int iCol = 0; iCol < outUpstreamR.size(); ++iCol){
              mask(outUpstreamR[iCol],outUpstreamR[iCol])=0;
            }
            downstream = mu*lR*mask;	
            
            // consider axial boundary values defined in this cell
            mask.setIdentity();
            for (int iCol = 0; iCol < outUpstreamZ.size(); ++iCol){
              mask(outUpstreamZ[iCol],outUpstreamZ[iCol])=0;
            }
            downstream = downstream + xi*lZ*mask;	
            
            A = A + downstream;

            // form b matrix
            b = t*q;
          
            // consider contribution of angular redistribution term 
            // calculated with known values
            angRedistCoeff = ((alphaPlusOneHalf/tau)*(tau - 1.0)\
              - alphaMinusOneHalf)/weight;
            cellHalfAFlux.setOnes();
            cellHalfAFlux=(*halfAFlux)(iZ,iR,iXi)*cellHalfAFlux;
            b = b - angRedistCoeff*R*cellHalfAFlux;

            // consider radial boundary values in other cells or BCs
            if (iR==rStart){
              // apply reflecting boundary condition

              // calculate indices of reflected ordinate
              numPs = mesh->quadrature.size();
              numQs = mesh->quadrature[iXi].nOrd;
              reflectedP= numPs-iXi-1;
              reflectedQ= numQs-iMu-1;

              // get angular index of reflect flux
              reflectedAngIdx = mesh->quadrature[reflectedP].ordIdx[reflectedQ]; 
              
              // account for reflected flux at boundary
              upstream = mu*(*aFlux)(iZ,iR,reflectedAngIdx)\
              *(lR.col(outUpstreamR[0])+lR.col(outUpstreamR[1]));
              b = b - upstream;
            }

            // consider radial boundary values defined in other cells
            // or by BCs
            else if (iR!=rStart){
              upstream = mu*(*aFlux)(iZ,iR+borderCellR,angIdx)\
              *(lR.col(outUpstreamR[0])+lR.col(outUpstreamR[1]));
              b = b - upstream;
            }

            // consider axial boundary values defined in other cells
            // or by BCs
            if (iZ!=zStart){
              upstream = xi*(*aFlux)(iZ+borderCellZ,iR,angIdx)\
              *(lZ.col(outUpstreamZ[0])+lZ.col(outUpstreamZ[1]));
              b = b - upstream;
            }

            // solve for angular fluxes in each corner 
            x = A.lu().solve(b);

            // take average of corner values to get angular flux
            // for this cell
            subCellVol = calcSubCellVol(iZ,iR);	
            (*aFlux)(iZ,iR,angIdx) = x.dot(subCellVol)/subCellVol.sum();

            // use weighted diamond difference to calculate next half
            // angle flux used for next nvalue of mu in this quadrature 
            // level
            (*halfAFlux)(iZ,iR,iXi) = ((*aFlux)(iZ,iR,angIdx)\
            +(tau-1.0)*(*halfAFlux)(iZ,iR,iXi))/tau; 

          } //iR
        } //iZ
      } // if mu > 0
    } //iMu
  } //iXi
  //cout << *aFlux << endl;
};
//==============================================================================

//==============================================================================
//! calckR calculate within cell radial leakage matrix

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
//! calckZ calculate within cell axial leakage matrix

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
//! calclR calculate out of cell radial leakage matrix

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
//! calclZ calculate out of cell axial leakage matrix

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
//! calct calculate first collision matrix

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
//! calcR calculate second collision matrix

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
//! calcSubCellVol calculate volumes of subcell regions

Eigen::VectorXd SimpleCornerBalance::calcSubCellVol(int myiZ, int myiR){
  Eigen::VectorXd subCellVol = Eigen::VectorXd::Zero(4);
          
  subCellVol(0) = (mesh->dzs(myiZ)/2)*(pow(mesh->rCent(myiR),2)-\
    pow(mesh->rEdge(myiR),2));

  subCellVol(1) = (mesh->dzs(myiZ)/2)*(pow(mesh->rEdge(myiR+1),2)-\
    pow(mesh->rCent(myiR),2));

  subCellVol(2) = (mesh->dzs(myiZ)/2)*(pow(mesh->rEdge(myiR+1),2)-\
    pow(mesh->rCent(myiR),2));

  subCellVol(3) = (mesh->dzs(myiZ)/2)*(pow(mesh->rCent(myiR),2)-\
    pow(mesh->rEdge(myiR),2));

  return subCellVol;
}
//==============================================================================
