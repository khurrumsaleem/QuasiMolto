// File: StartingAngle.cpp
// Purpose: calculate initial half angle fluxes needed for 
// approximation of the angular redistribution term of the
// neutron transport equation in RZ geometry
// Date: October 28, 2019

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Materials.h"
#include "Material.h"
#include "StartingAngle.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
/// StartingAngle object constructor
///
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myInput YAML input object for the simulation
StartingAngle::StartingAngle(Mesh * myMesh,\
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
/// Calculate the starting half angle angular flux
///
/// Solves using simple corner balance on a simplified neutron transport equation 
/// in RZ geometry. The simplification comes about by substituting quadrature 
/// values that correspond to the starting half angle, and results in the 
/// elimination of the angular redistribution term. 
void StartingAngle::calcStartingAngle(cube * halfAFlux,\
  Eigen::MatrixXd * source,\
  Eigen::MatrixXd * alpha,\
  int energyGroup)
{
  //cout << "Entered StartingAngle solve" << endl; 

  // Index xi value is stored in in quadLevel
  const int xiIndex = 0;

  // Temporary variable used for looping though quad set
  double xi,sqrtXi,sigTEff,sigTEps=1E-4;
  int zStart,rStart,zEnd,zInc,borderCellZ,borderCellR,zStartCell,rStartCell;
  int rows = 4,cols = 4;
  double v = materials->neutV(energyGroup);
  vector<int> withinUpstreamR(2);
  vector<int> outUpstreamR(2);
  vector<int> withinUpstreamZ(2);
  vector<int> outUpstreamZ(2);
  Eigen::MatrixXi cornerOffset(4,2);
  Eigen::MatrixXd sigT(rows,cols);
  double gamma;

  // Within cell leakage matrices in R and Z directions
  double kRCoeff;
  double kZCoeff;
  Eigen::MatrixXd kR = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd kZ = Eigen::MatrixXd::Zero(rows,cols);

  // Out of cell leakage matrices in R and Z directions
  double lRCoeff;
  double lZCoeff;
  Eigen::MatrixXd lR = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd lZ = Eigen::MatrixXd::Zero(rows,cols);

  // Reaction matrices
  double t1Coeff;
  double t2Coeff;
  Eigen::MatrixXd t1 = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd t2 = Eigen::MatrixXd::Zero(rows,cols);

  // A matrix of linear system Ax=b
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::MatrixXd mask = Eigen::MatrixXd::Zero(rows,cols);
  Eigen::VectorXd subCellVol = Eigen::VectorXd::Zero(rows);

  // Downstream values
  Eigen::MatrixXd downstream = Eigen::MatrixXd::Zero(rows,cols);

  // Downstream values
  Eigen::VectorXd upstream = Eigen::VectorXd::Zero(rows);

  // Right-hand side
  Eigen::VectorXd b = Eigen::VectorXd::Zero(rows);

  // Solution vector 
  Eigen::VectorXd x = Eigen::VectorXd::Zero(rows);

  // Source 
  Eigen::VectorXd q = Eigen::VectorXd::Ones(rows);

  // Dirichlet boundary condition
  double rBC,zBC;

  // MMS source 
  Eigen::VectorXd mmsQ = Eigen::VectorXd::Zero(rows);
  
  for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi){

    // Get xi for this quadrature level
    xi = mesh->quadrature[iXi].quad[0][xiIndex];
    sqrtXi = pow(1-pow(xi,2),0.5); 
    rStart = mesh->drsCorner.size()-1;
    rStartCell = mesh->drs.size()-1;
    borderCellR = 1;
    withinUpstreamR = {0,3};
    outUpstreamR = {1,2};
  
    // Set dirichlet bc
    rBC=outerBC[energyGroup];

    // Depending on xi, define loop constants	
    if (xi > 0) {
      // Marching from bottom to top
      zStart = 0;
      zStartCell = 0;
      zInc = 2;
      zEnd = 0;
      borderCellZ = -1;
      withinUpstreamZ = {2,3};
      outUpstreamZ = {0,1};
      cornerOffset << -1,0,\
        0,0,\
        0,1,\
        -1,1;
      

      // Set dirichlet bc
      zBC = lowerBC[energyGroup];
            
    }
    else {			
      
      // Marching from the top to the bottom  
      zStart = mesh->dzsCorner.size()-1;
      zStartCell = mesh->dzs.size()-1;
      zEnd = mesh->dzs.size();
      zInc = -2;
      borderCellZ = 1;
      withinUpstreamZ = {0,1};
      outUpstreamZ = {2,3};
      cornerOffset << -1,-1,\
        0,-1,\
        0,0,\
        -1,0;

      // Set dirichlet bc
      zBC = upperBC[energyGroup];
    }
    for (int iR = rStart,iCellR = rStartCell,countR = 0;\
      countR < mesh->drsCorner.size(); 
      iR = iR - 2,--iCellR,countR = countR + 2){

      for (int iZ = zStart, iCellZ = zStartCell,countZ = 0; \
        countZ < mesh->dzsCorner.size();\
        iZ = iZ + zInc, iCellZ = iCellZ + zInc/2,countZ = countZ + 2){
        
        // Set source in each corner
        for (int iCorner = 0; iCorner < 4; ++iCorner){
          q(iCorner) = (*source)(iZ+cornerOffset(iCorner,1),\
          iR+cornerOffset(iCorner,0));
        }

      //  sigT = materials->sigT(iZ,iR,energyGroup) + (*alpha)(iZ,iR)/v;
      //  if (sigT < sigTEps){
      //    sigT = sigTEps;
      //  }

        for (int iSig = 0; iSig < sigT.cols(); ++iSig){
          sigTEff = materials->sigT(iZ+cornerOffset(iSig,1),\
            iR+cornerOffset(iSig,0),energyGroup)\
            +(*alpha)(iZ+cornerOffset(iSig,1),iR+cornerOffset(iSig,0))/v;

          if (sigTEff > sigTEps)
            sigT(iSig,iSig) = sigTEff;
          else
            sigT(iSig,iSig) = sigTEps;
        }

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
        t1Coeff = mesh->drs(iCellR)*mesh->dzs(iCellZ)*mesh->rEdge(iCellR+1)/16.0;
        t1=calct1(gamma);
        t1=t1Coeff*t1;
        
        // Calculate second collision matrix
        t2Coeff = mesh->drs(iCellR)*mesh->dzs(iCellZ)/4.0;
        t2=calct2(gamma);
        t2=t2Coeff*t2;

        // Calculate corner volumes
        subCellVol = calcSubCellVol(iCellZ,iCellR);	

        // Calculate A considering within cell leakage and 
        // collision matrices
        A = sqrtXi*kR+xi*kZ+sigT*t1+sqrtXi*t2;

        // Consider radial downstream values defined in this cell
        mask.setIdentity();
        for (int iCol = 0; iCol < outUpstreamR.size(); ++iCol){
          mask(outUpstreamR[iCol],outUpstreamR[iCol])=0;
        }
        downstream = sqrtXi*lR*mask;	
        
        // Consider axial downstream values defined in this cell
        mask.setIdentity();
        for (int iCol = 0; iCol < outUpstreamZ.size(); ++iCol){
          mask(outUpstreamZ[iCol],outUpstreamZ[iCol])=0;
        }
        downstream = downstream + xi*lZ*mask;	
        
        A = A + downstream;
        
        // Form b matrix
        b = t1*q;
        // Consider upstream values in other cells or BCs
        if (iR!=rStart){
          //upstream = sqrtXi*(*halfAFlux)(iZ,iR+borderCellR,iXi)\
          *(lR.col(outUpstreamR[0])+lR.col(outUpstreamR[1]));
          upstream =\
          sqrtXi*(*halfAFlux)(iZ+cornerOffset(outUpstreamR[0],1),\
          iR+borderCellR,iXi)*lR.col(outUpstreamR[0])+\
          sqrtXi*(*halfAFlux)(iZ+cornerOffset(outUpstreamR[1],1),\
          iR+borderCellR,iXi)*lR.col(outUpstreamR[1]);
          b = b - upstream;
        } else {
          upstream = sqrtXi*rBC\
          *(lR.col(outUpstreamR[0])+lR.col(outUpstreamR[1]));
          b = b - upstream;
        }
        if (iZ!=zStart){
          upstream =\
          xi*(*halfAFlux)(iZ+borderCellZ,iR+cornerOffset(outUpstreamZ[0],0),\
          iXi)*lZ.col(outUpstreamZ[0])+\
          xi*(*halfAFlux)(iZ+borderCellZ,iR+cornerOffset(outUpstreamZ[1],0),\
          iXi)*lZ.col(outUpstreamZ[1]);
          b = b - upstream;
        }else{
          upstream = xi*zBC\
          *(lZ.col(outUpstreamZ[0])+lZ.col(outUpstreamZ[1]));
          b = b - upstream;
        }
       
        mmsQ = calcMMSSource(iCellZ,iCellR,energyGroup,iXi,sigT,subCellVol);
        b = b + 0.25*t1*mmsQ;

        x = A.partialPivLu().solve(b);
        
        // Take average of half angle fluxes in corners
        //(*halfAFlux)(iZ,iR,iXi) = x.dot(subCellVol)/subCellVol.sum();
        //(*halfAFlux)(iZ,iR,iXi) = x.sum()/4.0;
        
        (*halfAFlux)(iZ+cornerOffset(0,1),iR+cornerOffset(0,0),iXi) = x(0);
        (*halfAFlux)(iZ+cornerOffset(1,1),iR+cornerOffset(1,0),iXi) = x(1);
        (*halfAFlux)(iZ+cornerOffset(2,1),iR+cornerOffset(2,0),iXi) = x(2);
        (*halfAFlux)(iZ+cornerOffset(3,1),iR+cornerOffset(3,0),iXi) = x(3);
      }
    }
  }
  cout << "half angle angular flux" << endl;
  cout << (*halfAFlux) << endl;	
};
//==============================================================================

//==============================================================================
/// Calculate within cell radial leakage matrix
///
/// @param [in] myGamma Ratio of inner radius to outer radius in a cell
/// @param [out] kR Within cell radial leakage matrix 
Eigen::MatrixXd StartingAngle::calckR(double myGamma){
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
/// Calculate within cell axial leakage matrix
///
/// @param [in] myGamma Ratio of inner radius to outer radius in a cell
/// @param [out] kZ Within cell axial leakage matrix 
Eigen::MatrixXd StartingAngle::calckZ(double myGamma){
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
Eigen::MatrixXd StartingAngle::calclR(double myGamma){
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
/// Calculate out of cell axial leakage matrix
///
/// @param [in] myGamma Ratio of inner radius to outer radius in a cell
/// @param [out] lZ Out of cell axial leakage matrix 
Eigen::MatrixXd StartingAngle::calclZ(double myGamma){
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
/// Calculate first collision matrix
///
/// @param [in] myGamma Ratio of inner radius to outer radius in a cell
/// @param [out] t1 First collision matrix 
Eigen::MatrixXd StartingAngle::calct1(double myGamma){
  double a = 1+3*myGamma;
  double b = 3+myGamma;
  Eigen::MatrixXd t1 = Eigen::MatrixXd::Zero(4,4);

  t1(0,0) = a; 
  t1(1,1) = b; 
  t1(2,2) = b; 
  t1(3,3) = a; 
   
  return t1;
}
//==============================================================================

//==============================================================================
/// Calculate second collision matrix
///
/// @param [in] myGamma Ratio of inner radius to outer radius in a cell
/// @param [out] t2 Second collision matrix 
Eigen::MatrixXd StartingAngle::calct2(double myGamma){
  double a = 1;
  Eigen::MatrixXd t2 = Eigen::MatrixXd::Zero(4,4);

  t2(0,0) = a; 
  t2(1,1) = a; 
  t2(2,2) = a;
  t2(3,3) = a; 
   
  return t2;
}
//==============================================================================

//==============================================================================
/// Calculate volumes of subcell regions
///
/// @param [out] subCellVol Volume in each corner of a cell
Eigen::VectorXd StartingAngle::calcSubCellVol(int myiZ, int myiR){
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

//==============================================================================
/// Calculate MMS source in subcell regions
///
/// @param [out]  mmsSource Contains the mmsSource in each cell
Eigen::VectorXd StartingAngle::calcMMSSource(int myiZ,int myiR,\
  int energyGroup, int iXi, Eigen::MatrixXd sigT, Eigen::VectorXd subCellVol){
 
  // Define material parameters
  double sigS = materials->sigS(myiZ,myiR,energyGroup,energyGroup);
  double sigF = materials->sigF(myiZ,myiR,energyGroup);
  double mySigT = materials->sigT(myiZ,myiR,energyGroup);
  double nu = materials->nu(myiZ,myiR);
  double xi = mesh->quadrature[iXi].quad[0][0]; 
  double mu = -1.0*sin(acos(xi)); 
  double eta = 0.0; 
  double t = 0.0001;
  double gamma = sin(acos(xi));
 
  // Return variable
  Eigen::VectorXd mmsSource = Eigen::VectorXd::Zero(4);
  
  // Define bounds corners are integrated over
  Eigen::MatrixXd z = Eigen::MatrixXd::Zero(2,2);
  z << mesh->zCent(myiZ), mesh->zEdge(myiZ+1),\
    mesh->zEdge(myiZ),mesh->zCent(myiZ); 
  z << mesh->zEdge(myiZ),mesh->zCent(myiZ),
    mesh->zCent(myiZ), mesh->zEdge(myiZ+1);
  Eigen::MatrixXd r = Eigen::MatrixXd::Zero(4,2);
  r << mesh->rEdge(myiR),mesh->rCent(myiR),\
    mesh->rCent(myiR), mesh->rEdge(myiR+1),\
    mesh->rCent(myiR), mesh->rEdge(myiR+1),\
    mesh->rEdge(myiR),mesh->rCent(myiR);
 // z << mesh->zEdge(myiZ),mesh->zEdge(myiZ+1),\
    mesh->zEdge(myiZ), mesh->zEdge(myiZ+1);
//  r << mesh->rEdge(myiR),mesh->rEdge(myiR+1),\
    mesh->rEdge(myiR), mesh->rEdge(myiR+1);
  // Get maximum radius and axial height
  double R = mesh->R,Z = mesh->Z;
  
  // Set coefficient in exponential
  double c = 1.0;
  double v = 1.0;

  // Define some temporary variables to clean things up
  double A = exp(c*t) * pow(mu,3);
  double B = exp(c*t) * M_PI * xi * pow(mu,2) / Z;
  double D = -exp(c*t) * mu * ( pow(sin(acos(xi)),2) - 3*pow(eta,2));
  double F = exp(c*t) * ( pow(mu,2) * (mySigT + c/v) - (4.0*M_PI/3.0)\
    *(sigS + nu * sigF)); 

  A = 2.0*exp(c*t) * pow(sin(gamma),3);
  B = exp(c*t) * M_PI * xi * pow(sin(gamma),2) / Z;
  D = exp(c*t) * ( pow(sin(gamma),2) * (mySigT + c/v) - (4.0*M_PI/3.0)\
    *(sigS + nu * sigF)); 

  A = exp(c*t) * (mySigT + c/v -1.0*(sigS + nu * sigF));
  B = exp(c*t) * mu;
  D = exp(c*t) * M_PI * xi/ Z;
 
  // Temporary variables to evaluate whole MMS term
  double term1,term2,term3,term4,term5;
  double rad1,rad2,rad3,rad4;
  double sinUp,sinDown,cosUp,cosDown;

  // Counter to index return variable 
  int count = 0; 

  // Loop over each corner and calculate source
  for (int iAx = 0; iAx < 2; ++iAx){
    for (int iRad = 0; iRad < 2; ++iRad){
    
      // Evaluate sine and cosine terms
      sinUp = sin(M_PI*z(iAx,1)/Z);
      sinDown = sin(M_PI*z(iAx,0)/Z);
      cosUp = cos(M_PI*z(iAx,1)/Z);
      cosDown = cos(M_PI*z(iAx,0)/Z);
      
      // Evaluate constituent terms
      rad1 = (pow(r(iRad,1),3)/3.0)- (pow(r(iRad,0),3)/3.0);
      term1 = -(A*Z/M_PI)*(cosUp-cosDown)*rad1;

      rad2 = ((pow(R,2)*pow(r(iRad,1),2)/2.0)-(pow(r(iRad,1),4)/4.0))\
        - ((pow(R,2)*pow(r(iRad,0),2)/2.0)-(pow(r(iRad,0),4)/4.0));
      term2 = (B*Z/M_PI)*(sinUp-sinDown)*rad2;

      rad3=rad2;
      term3 = -(D*Z/M_PI)*(cosUp-cosDown)*rad3;


      // Evaluate constituent terms
 //     rad1 = ((pow(R,2)*r(iRad,1)-pow(r(iRad,1),3))\
  //      -(pow(R,2)*r(iRad,0) - pow(r(iRad,0),3)));
  //    term1 = -(A*Z/M_PI)*(cosUp-cosDown)*rad1;

    //  rad2 = ((pow(R,2)*pow(r(iRad,1),2)/2.0-pow(r(iRad,1),4)/4.0)
    //    -(pow(R,2)*pow(r(iRad,0),2)/2.0 - pow(r(iRad,0),4)/4.0));
    //  term2 = (B*Z/M_PI)*(sinUp-sinDown)*rad2;

     // rad3 = ((pow(R,2)*r(iRad,1)-pow(r(iRad,1),3)/3.0)\
     //   -(pow(R,2)*r(iRad,0)-pow(r(iRad,0),3)/3.0));
     // term3 = -(D*Z/M_PI)*(cosUp-cosDown)*rad3;

     // rad4 = ((pow(R,2)*pow(r(iRad,1),2)/2.0-pow(r(iRad,1),4)/4.0)\
     //   -(pow(R,2)*pow(r(iRad,0),2)/2.0-pow(r(iRad,0),4)/4.0));
    //  term4 = -(F*Z/M_PI)*(cosUp-cosDown)*rad4;

      // define angular independent mms
     // mmsSource(count) = -8.0 * M_PI * (term1 - term2 + term3 + term4);
    //  mmsSource(count) = 4.0 * (term1 + term2 + term3 + term4);
    //  ++count;
      rad1 = ((pow(R,2)*pow(r(count,1),2)/2.0-pow(r(count,1),4)/4.0)
        -(pow(R,2)*pow(r(count,0),2)/2.0 - pow(r(count,0),4)/4.0));
      term1 = -(A*Z/M_PI)*(cosUp-cosDown)*rad1;

      rad2 = ((pow(R,2)*r(count,1)-pow(r(count,1),3)/3.0)\
        -(pow(R,2)*r(count,0)-pow(r(count,0),3)/3.0));
      term2 = (B*Z/M_PI)*(cosUp-cosDown)*rad2;

      rad3 = ((pow(R,2)*r(count,1)-pow(r(count,1),3))\
        -(pow(R,2)*r(count,0) - pow(r(count,0),3)));
      term3 = -(B*Z/M_PI)*(cosUp-cosDown)*rad3;

      rad4 = ((pow(R,2)*pow(r(count,1),2)/2.0-pow(r(count,1),4)/4.0)\
        -(pow(R,2)*pow(r(count,0),2)/2.0-pow(r(count,0),4)/4.0));
      term4 = (D*Z/M_PI)*(sinUp-sinDown)*rad4;

      term5 = mySigT*subCellVol(count);
      term5 = 0.0;

      mmsSource(count) = 4.0* (term1 + term2 + term3 + term4 + term5)/subCellVol(count);
//      mmsSource(count) = 1.0/subCellVol(count); 
//      mmsSource(count) = 4.0* (term1 + term2 + term3 + term4) * \
        subCellVol(count)/subCellVol.sum() + 4.0*term5;
      ++count;
    }
  }
  
        
  return mmsSource;
}
//==============================================================================


