// File: StartingAngle.cpp
// Purpose: calculate initial half angle fluxes needed for 
// approximation of the angular redistribution term of the
// neutron transport equation in RZ geometry
// Date: October 28, 2019

#include "StartingAngle.h"

using namespace std; 

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

  // ToDo: make function to read in optional arguments
  // Check for optional inputs
  if ((*input)["parameters"]["upperAngularFluxBC"]){
    inpUpperBC=(*input)["parameters"]["upperAngularFluxBC"].as<vector<double>>();

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

  if ((*input)["parameters"]["lowerAngularFluxBC"]){
    inpLowerBC=(*input)["parameters"]["lowerAngularFluxBC"].as<vector<double>>();

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

  if ((*input)["parameters"]["outerAngularFluxBC"]){
    inpOuterBC=(*input)["parameters"]["outerAngularFluxBC"].as<vector<double>>();

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
void StartingAngle::calcStartingAngle(arma::cube * halfAFlux,\
    Eigen::MatrixXd * source,\
    Eigen::MatrixXd * alpha,\
    int energyGroup)
{

  for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi){
    solveAngularFlux(halfAFlux,source,alpha,energyGroup,iXi);
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
void StartingAngle::solveAngularFlux(arma::cube * halfAFlux,\
    Eigen::MatrixXd * source,\
    Eigen::MatrixXd * alpha,\
    int energyGroup,int iXi)
{


  // Index xi value is stored in in quadLevel
  const int xiIndex = 0;

  // Temporary variable used for looping though quad set
  double xi,sqrtXi,sigTEff,v,sigTEps=1E-4;

  xi = mesh->quadrature[iXi].quad[0][xiIndex];
  int zStart,rStart,zEnd,zInc,borderCellZ,borderCellR,zStartCell,rStartCell;
  int rows = 4,cols = 4;
  vector<int> withinUpstreamR(2);
  vector<int> outUpstreamR(2);
  vector<int> withinUpstreamZ(2);
  vector<int> outUpstreamZ(2);
  Eigen::MatrixXi cornerOffset(4,2);
  Eigen::MatrixXd sigT = Eigen::MatrixXd::Zero(rows,cols);
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
  Eigen::VectorXd q = Eigen::VectorXd::Zero(rows);

  // Dirichlet boundary condition
  double rBC,zBC;

  // MMS source 
  Eigen::VectorXd mmsQ = Eigen::VectorXd::Zero(rows);

  // Get xi for this quadrature level
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

      for (int iSig = 0; iSig < sigT.cols(); ++iSig){
        
        // Get neutron velocity in this corner
        v = materials->neutVel(iZ,iR,energyGroup);

        // Calculate effective cross section in this corner 
        v = materials->neutVel(iZ,iR,energyGroup);
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

      x = A.partialPivLu().solve(b);

      (*halfAFlux)(iZ+cornerOffset(0,1),iR+cornerOffset(0,0),iXi) = x(0);
      (*halfAFlux)(iZ+cornerOffset(1,1),iR+cornerOffset(1,0),iXi) = x(1);
      (*halfAFlux)(iZ+cornerOffset(2,1),iR+cornerOffset(2,0),iXi) = x(2);
      (*halfAFlux)(iZ+cornerOffset(3,1),iR+cornerOffset(3,0),iXi) = x(3);
    }
  }
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


