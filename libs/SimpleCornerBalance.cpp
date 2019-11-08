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
//! solve Calculate the starting half angle angular flux

// Solves using simple corner balance on the neutron transport equation 
// in RZ geometry. 

void SimpleCornerBalance::solve(cube * aFlux,\
  cube * halfAFlux,\
  Eigen::MatrixXd * source,\
  int energyGroup)
{
  // index xi value is stored in in quadLevel
  const int xiIndex = 0;
  // index xi value is stored in in quadLevel
  const int muIndex = 1;
  // index weight value is stored in in quadLevel
  const int weightIndex = 3;
  // temporary variable used for looping though quad set
  double xi,sqrtXi,sigT,mu,alphaMinusOneHalf,alphaPlusOneHalf,weight,tau;
  double angRedistCoeff = 0;
  int zStart,rStart,zEnd,zInc,borderCellZ,borderCellR,angIdx;
  int rows = 4,cols = 4;
  vector<int> withinUpstreamR(2);
  vector<int> outUpstreamR(2);
  vector<int> withinUpstreamZ(2);
  vector<int> outUpstreamZ(2);
  double gamma;
  // within cell leakage matrices in R and Z directions
  double kRCoeff;
  double kZCoeff;
  Eigen::MatrixXd kR = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd kZ = Eigen::MatrixXd::Zero(rows,cols);
  // out of cell leakage matrices in Rand Z directions
  double lRCoeff;
  double lZCoeff;
  Eigen::MatrixXd lR = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd lZ = Eigen::MatrixXd::Zero(rows,cols);
  // reaction matrices
  double tCoeff;
  double RCoeff;
  Eigen::MatrixXd t = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(rows,cols);
  // A matrix of linear system Ax=b
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd mask = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::VectorXd subCellVol = Eigen::VectorXd::Zero(rows);
  // downstream values
  Eigen::MatrixXd downstream = Eigen::MatrixXd::Zero(rows,cols);
  // downstream values
  Eigen::VectorXd upstream = Eigen::VectorXd::Zero(rows);
  // half-angle fluxes 
  Eigen::VectorXd cellHalfAFlux = Eigen::VectorXd::Zero(rows);
  // right-hand side
  Eigen::VectorXd b = Eigen::VectorXd::Zero(rows);
  // solution vector 
  Eigen::VectorXd x = Eigen::VectorXd::Zero(rows);
  // source 
  Eigen::VectorXd q = Eigen::VectorXd::Ones(rows);
  q = 8*q;

  for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi){
    // get xi for this quadrature level
    xi = mesh->quadrature[iXi].quad[0][xiIndex];

    for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu){

      // get xi for this ordinate 
      mu = mesh->quadrature[iXi].quad[iMu][muIndex]; 
      alphaPlusOneHalf = mesh->quadrature[iXi].alpha[iMu+1];
      alphaMinusOneHalf = mesh->quadrature[iXi].alpha[iMu];
      weight = mesh->quadrature[iXi].quad[iMu][weightIndex]; 
      tau = mesh->quadrature[iXi].tau[iMu];
      angIdx= mesh->quadrature[iXi].ordIdx[iMu];

      if (mu < 0 ){

        rStart = mesh->drs.size()-1;
        borderCellR = 1;
        withinUpstreamR = {0,3};
        outUpstreamR = {1,2};

        // depending on xi, define loop constants	
        if (xi > 0) {
          zStart = mesh->dzs.size()-1;
          zEnd = 0;
          zInc = -1;
          borderCellZ = 1;
          withinUpstreamZ = {2,3};
          outUpstreamZ = {0,1};
        }
        else {			
          zStart = 0;
          zEnd = mesh->dzs.size();
          zInc = 1;
          borderCellZ = -1;
          withinUpstreamZ = {0,1};
          outUpstreamZ = {2,3};
        }
        for (int iR = rStart, countR = 0; countR < mesh->drs.size(); --iR, ++countR){
          for (int iZ = zStart, countZ = 0; \
            countZ < mesh->dzs.size(); iZ = iZ + zInc, ++countZ){

            q.setOnes();
            q = (*source)(iZ,iR)*q;
            sigT = materials->sigT(iZ,iR,energyGroup);
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

            // calculate A considering within cell leakage and 
            // collision matrices
            angRedistCoeff = alphaPlusOneHalf/(weight*tau); 
            A = mu*kR+xi*kZ+sigT*t+angRedistCoeff*R;

            // consider radial downstream values defined in this cell
            mask.setIdentity();
            for (int iCol = 0; iCol < outUpstreamR.size(); ++iCol){
              mask(outUpstreamR[iCol],outUpstreamR[iCol])=0;
            }
            downstream = mu*lR*mask;	
            
            // consider axial downstream values defined in this cell
            mask.setIdentity();
            for (int iCol = 0; iCol < outUpstreamZ.size(); ++iCol){
              mask(outUpstreamZ[iCol],outUpstreamZ[iCol])=0;
            }
            downstream = downstream + xi*lZ*mask;	
            
            A = A + downstream;

            // form b matrix
            b = t*q;

            angRedistCoeff = ((alphaPlusOneHalf/tau)*(tau - 1.0)\
              - alphaMinusOneHalf)/weight;
            cellHalfAFlux.setOnes();
            cellHalfAFlux=(*halfAFlux)(iZ,iR,iXi)*cellHalfAFlux;
            b = b - angRedistCoeff*R*cellHalfAFlux;

            // consider upstream values in other cells or BCs
            if (iR!=rStart){
              upstream = mu*(*halfAFlux)(iZ,iR+borderCellR,iXi)\
              *(lR.col(outUpstreamR[0])+lR.col(outUpstreamR[1]));
              b = b - upstream;
            }

            if (iZ!=zStart){
              upstream = xi*(*halfAFlux)(iZ+borderCellZ,iR,iXi)\
              *(lZ.col(outUpstreamZ[0])+lZ.col(outUpstreamZ[1]));
              b = b - upstream;
            }

            x = A.lu().solve(b);

            // Take average of subcells
            subCellVol = calcSubCellVol(iZ,iR);	
            (*aFlux)(iZ,iR,angIdx) = x.dot(subCellVol)/subCellVol.sum();

            // use weighted diamond difference to calculate next half
            // angle flux
            (*halfAFlux)(iZ,iR,iXi) = ((*aFlux)(iZ,iR,angIdx)\
            +(tau-1.0)*(*halfAFlux)(iZ,iR,iXi))/tau; 

          } //iR
        } //iZ
      } else if (mu > 0){
        //case when mu > 0
      }
    } //iMu
  } //iXi
	
  cout << "angular flux calculated! " << endl;
  cout << *aFlux << endl;
};
//==============================================================================

//==============================================================================
//! calckR calculate within cell radial leakage matrix

Eigen::MatrixXd SimpleCornerBalance::calckR(double myGamma){
  double a = -(1+myGamma);
  double b = 1+myGamma;
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
  double a = myGamma;
  double b = -1;
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
