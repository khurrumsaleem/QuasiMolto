// File: mesh.cpp
// Purpose: defines a class to hold mesh data such as spatial mesh and 
// quadrature set
// Author: Aaron James Reynolds
// Date: October 9, 2019

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"

using namespace std;
using namespace arma;

//==============================================================================
/// Constructor for quadLevel object.
///
/// @param [in] myQuad Contains ordinates on a single quadrature level 
/// @param [in] myAlpha Differencing coefficients for a single quadrature level 
/// @param [in] myTau Coefficient for weighted diamond differencing on a single 
/// quadrature level
/// @param [in] myStartIndex Starting index for ordinates on the quadrature level
/// to be created 
quadLevel::quadLevel(vector< vector<double> > myQuad,\
  vector<double> myAlpha,\
  vector<double> myTau,\
  int myStartIndex){

  int ordIdxCount = myStartIndex;	

  // Resize vectors on quadLevel
  quad.resize(myQuad.size(),vector<double>(myQuad[0].size(),0.0));
  alpha.resize(myAlpha.size());
  tau.resize(myTau.size());
  ordIdx.resize(myQuad.size());
  
  // Set equal to initializing arguments
  quad = myQuad;
  alpha = myAlpha;
  tau = myTau;
  nOrd = quad.size();
  
  // Assign index to each ordinate
  for (int iOrd = 0; iOrd < nOrd; ++iOrd){
    ordIdx[iOrd]=ordIdxCount;
    ++ordIdxCount;
  }
}
//==============================================================================

//==============================================================================
/// Mesh Contructor for Mesh object.
///
/// @param [in] myInput YAML input object for this simulation 
Mesh::Mesh(YAML::Node * myInput){

  // Set input pointer
  input = myInput;

  // Currently only works for a quadrature set of order 12 
  n=12; 

  // Read in input mesh parameters
  dz = (*input)["mesh"]["dz"].as<double>();
  dr = (*input)["mesh"]["dr"].as<double>();
  Z = (*input)["mesh"]["Z"].as<double>();
  R = (*input)["mesh"]["R"].as<double>();
  dt = (*input)["mesh"]["dt"].as<double>();
  
  // Set up the mesh and quadrature set
  calcSpatialMesh();
  calcQuadSet();
  calcNumAnglesTotalWeight();

}
//==============================================================================


//==============================================================================
/// Calculates a quadrature set and differencing coefficients
///
/// Uses the methodology laid out Methods of Computational Transport by
/// Lewis and Miller. A level symmetric quadrature set is calculated. The p,q 
/// indexing scheme is defined on page 166 in Figure 4-7. Differencing 
/// coefficients are then calculated using that quadrature set.
/// TODO: a lot of this is very specific to an n=12 quadrature set. Would be 
/// better if it was more general          
void Mesh::calcQuadSet(){	
	
  // Weights in quadrature set given in L&M
  double weight1 = 0.0707626;
  double weight2 = 0.0558811;
  double weight3 = 0.0373377;
  double weight4 = 0.0502819;
  double weight5 = 0.0258513;

  // Positions corresponding to each weight in L&M
  vector<int> ord1 = {1,6,21}; 
  vector<int> ord2 = {2,5,7,11,19,20};
  vector<int> ord3 = {3,4,12,15,16,18};
  vector<int> ord4 = {8,10,17};
  vector<int> ord5 = {9,13,14};

  // temporary variable for building ascending quadrature list
  vector< vector<double> > tempQuad;
  vector< vector<double> > orderedQuad;
  vector< vector<double> > alpha; // differencing coefficients
  vector<double> tempRow;
  vector<int> coeffQ1;
  vector<int> coeffQ2;
  int counter,rowIndexOfMinimum;

  // calculate viable mu according to L&M approach
  calcMu();

  // Find number of sets of ordinates whose magnitude is one in the first
  // quadrant
  counter = 0;
  for (int i = 0; i < n/2; ++i){
    for (int j =0; j < n/2; ++j){
      for (int k =0; k < n/2; ++k){
        if (abs(pow(mu[i],2.0)+pow(mu[j],2.0)+pow(mu[k],2.0)-1)<1E-5){
          ++counter;
        }
      }
    }
  }

  // Size ordinates array with counter
  ordinates.resize(counter,vector<double>(4,0.0)); 
  counter = 0;

  // Fill ordinates with viable sets
  for (int i = 0; i < n/2; ++i){
    for (int j = 0; j < n/2; ++j){
      for (int k = 0; k < n/2; ++k){
        if (abs(pow(mu[i],2.0)+pow(mu[j],2.0)+pow(mu[k],2.0)-1)<1E-5){

        ordinates[counter][0] = mu[i];
        ordinates[counter][1] = mu[j];
        ordinates[counter][2] = mu[k];
        ++counter;
        
        }
      }
    }
  }

  // Assign each ordinate set the appropriate weight.
  for (int i = 0; i < ord1.size(); ++i)
    ordinates[ord1[i]-1][3] = weight1;

  for (int i = 0; i < ord2.size(); ++i)
    ordinates[ord2[i]-1][3] = weight2;

  for (int i = 0; i < ord3.size(); ++i)
    ordinates[ord3[i]-1][3] = weight3;

  for (int i = 0; i < ord4.size(); ++i)
    ordinates[ord4[i]-1][3] = weight4;

  for (int i = 0; i < ord5.size(); ++i)
    ordinates[ord5[i]-1][3] = weight5;

  // Now build quadrature over four quadrants
  coeffQ1 = {1,-1,1,-1};
  coeffQ2 = {1,1,-1,-1};
  tempQuad.resize(4*counter,vector<double>(4,0.0)); 
  for (int iquad = 0; iquad < 4; ++iquad){

    for (int iset = 0; iset < ordinates.size(); ++iset){

      tempRow = {coeffQ1[iquad]*ordinates[iset][0],\
        coeffQ2[iquad]*ordinates[iset][1],\
        ordinates[iset][2],ordinates[iset][3]};

      tempQuad[counter*iquad+iset] = tempRow;
    }
  }

  // Sort the quadrature set in ascending order	
  orderedQuad.resize(4*counter,vector<double>(3,0.0)); 

  for (int iOrdered = 0; iOrdered < orderedQuad.size(); ++iOrdered){
    rowIndexOfMinimum = 0;

    for (int iTemp = 0; iTemp < tempQuad.size(); ++iTemp){

      if (tempQuad[iTemp][0]<tempQuad[rowIndexOfMinimum][0]){
        rowIndexOfMinimum = iTemp;

      } else if (tempQuad[iTemp][0]==tempQuad[rowIndexOfMinimum][0]){

        if (tempQuad[iTemp][1]<tempQuad[rowIndexOfMinimum][1]){
          rowIndexOfMinimum = iTemp;	

        } else if (tempQuad[iTemp][1]==tempQuad[rowIndexOfMinimum][1]){

          if (tempQuad[iTemp][2]<tempQuad[rowIndexOfMinimum][2]){
            rowIndexOfMinimum = iTemp;
          }
        }
      }
    }

    // Place the next smallest on orderedQuad
    orderedQuad[iOrdered] = tempQuad[rowIndexOfMinimum];

    // Set this last placed row to some large values so we don't
    // accidentally recognize it as a minimum again
    tempQuad[rowIndexOfMinimum] = {100,100,100,100};
  }

  // Set quadSet
  quadSet.resize(4*counter,vector<double>(4,0.0)); 
  quadSet = orderedQuad;

  // Calculate differencing coefficients knowing the quadSet
  calcAlpha();
  calcTau();
  addLevels();
        
}
//==============================================================================

//==============================================================================
/// Calculate ordinates for a level symmetric quadrature set

void Mesh::calcMu(){

  // Allocate memory and fix only degree of freedom
  mu.resize(n/2,0.0);
  mu[0] = 0.1672126;

  // myConstant defined as on page 160 of L&M.
  double myConstant = 2.0*(1.0-3.0*pow(mu[0],2.0))/(n-2.0);
  for (int imu = 1; imu < n/2; ++imu){
          mu[imu] = sqrt(pow(mu[imu - 1],2.0) + myConstant);
  }
}
//==============================================================================

//==============================================================================
/// Calculate differencing coefficients for approximating angular redist term
///
/// Based on approach in Lewis and Miller
void Mesh::calcAlpha(){
  
  vector<int> rowLength = {2,4,6,8,10,12,12,10,8,6,4,2};
  alpha.resize(12,vector<double>(13,0.0)); 
  
  for (int i = 0; i < alpha.size(); ++i){
    
    // Initialize to 0 on the edge case to conserve neutrons. Explanation 
    // on page 179 of L&M.
    alpha[i][0] = 0; 
    
    // Use recursive definition to define subsequent values of alpha
    for (int j = 0; j < rowLength[i]; ++j){
      
      alpha[i][j+1] = alpha[i][j]-quadSet[quad_index(i,j)][1]\
        *quadSet[quad_index(i,j)][3];

      // Force sufficiently small quantities to 0
      // keeps negative values resulting from finite 
      // precision from coming up, too
      if (alpha[i][j+1]<1E-15)
        alpha[i][j+1] = 0.0;
    }
  }
}
//==============================================================================

//==============================================================================
/// Calculate tau value used in determining half angle angular fluxes
///
/// Based on approach in Lewis and Miller
void Mesh::calcTau(){
	
  double levelWeight = 0.0;
  vector<int> rowLength = {2,4,6,8,10,12,12,10,8,6,4,2};
  vector<double> halfOmega;
  vector<double> halfMu;

  tau.resize(12,vector<double>(12,0.0)); 
  for (int i = 0; i < tau.size(); ++i){
    
    // Calculate total weight on this level
    levelWeight = 0.0;
    for (int iWeight = 0; iWeight < rowLength[i]; ++iWeight){
      levelWeight=levelWeight+quadSet[quad_index(i,iWeight)][3];
    } 	
    
    // Calculate half angles on this level
    halfOmega.resize(rowLength[i]+1,0.0);
    halfOmega[0] = 0.0;

    for (int iHalfOmega = 1; iHalfOmega < halfOmega.size(); ++iHalfOmega){

      halfOmega[iHalfOmega] = halfOmega[iHalfOmega-1]\
      + M_PI*quadSet[quad_index(i,iHalfOmega-1)][3]/levelWeight;

    }

    // Calculate mu ordinates that correspond to each half angle
    halfMu.resize(rowLength[i]+1,0.0);
    halfMu[0] = sqrt(1-pow(quadSet[quad_index(i,0)][0],2))\
      *cos(halfOmega[0]);

    for (int iHalfMu = 1; iHalfMu < halfMu.size(); ++iHalfMu){

      halfMu[iHalfMu] = \
      sqrt(1-pow(quadSet[quad_index(i,iHalfMu-1)][0],2))\
      *cos(halfOmega[iHalfMu]);

    }
    
    // Calculate tau on this level
    for (int iTau = 0; iTau < rowLength[i]; ++iTau){
         
      tau[i][rowLength[i]-iTau-1] = \
      (quadSet[quad_index(i,rowLength[i]-1-iTau)][1]-halfMu[iTau])\
      /(halfMu[iTau+1]-halfMu[iTau]);
      
    }
  }
}
//==============================================================================

//==============================================================================
/// Build uniform mesh

void Mesh::calcSpatialMesh(){

  int nCellsZ;
  int nCellsR;

  // Calculate number of cells
  nCellsZ= Z/dz;
  nCellsR = R/dr;

  // Resize vectors storing dimensions of each cell
  dzs.zeros(2*nCellsZ);
  dzs.fill(dz/2.0);
  drs.zeros(2*nCellsR);
  drs.fill(dr/2.0);

  // Resize vector holding boundaries in each dimension
  rEdge.zeros(nCellsR+1);
  zEdge.zeros(nCellsZ+1);
  
  // Resize cell center location in each dimension
  rCent.zeros(nCellsR);
  zCent.zeros(nCellsZ);

  // Populate vectors holding boundaries in each dimension
  for (int iEdge = 1; iEdge < rEdge.n_elem; ++iEdge){
    rEdge(iEdge) = rEdge(iEdge-1) + drs(iEdge-1);
  }
  
  for (int iEdge = 1; iEdge < zEdge.n_elem; ++iEdge){
    zEdge(iEdge) = zEdge(iEdge-1) + dzs(iEdge-1);
  }
  
  // Populate vectors holding cell center location in each dimension
  for (int iCent = 0; iCent < rCent.n_elem; ++iCent){
    rCent(iCent) = (rEdge(iCent)+rEdge(iCent+1))/2.0;
  }
  
  for (int iCent = 0; iCent < zCent.n_elem; ++iCent){
    zCent(iCent) = (zEdge(iCent)+zEdge(iCent+1))/2.0;
  }

  // Calculate cell volume
  cellVol.resize(nCellsZ,vector<double>(nCellsR,0.0));

  for (int iZ = 0; iZ < cellVol.size(); ++iZ){
    for (int iR = 0; iR < cellVol[iZ].size(); ++iR){

      cellVol[iZ][iR] = dzs[iZ]*M_PI*\
      (pow(rEdge[iR+1],2)-pow(rEdge[iR],2));

    }
  }

  // Calculate area of vertical surfaces
  cellVSA.resize(nCellsZ,vector<double>(nCellsR+1,0.0));

  for (int iZ = 0; iZ < cellVSA.size(); ++iZ){
    for (int iR = 0; iR < cellVSA[iZ].size(); ++iR){

      cellVSA[iZ][iR] = dzs[iZ]*2*M_PI*rEdge[iR];

    }
  }

}
//==============================================================================

//==============================================================================
/// Add quadrature levels to the mesh object

void Mesh::addLevels(){

  vector< vector<double> > tempQuad;
  vector<double> tempAlpha;
  vector<double> tempTau;
  int count = 0;
  int startIndex = 0;
  int stopIndex = 0;
  int levelCount = 0;
  double myXi = 0.0;
  
  // Initialize first myXi to search for
  myXi=quadSet[0][0];
  for (int i = 0; i < quadSet.size(); ++i){

    // If the value of xi changes, set the stop index and copy 
    // the entries between start and stop indices
    if (myXi != quadSet[i][0]){
      stopIndex = i - 1;

      // Resize temporary variables
      tempQuad.resize(stopIndex-startIndex+1,vector<double>(4,0.0));
      tempAlpha.resize(stopIndex-startIndex+2,0.0);
      tempTau.resize(stopIndex-startIndex+1,0);
      
      // Counter for temporary arrays
      count = 0; 
      
      // Initialize boundary value of alpha
      tempAlpha[0] = alpha[levelCount][0];
      for (int iLevel = startIndex; iLevel < stopIndex + 1; ++iLevel){
        tempQuad[count] = quadSet[iLevel];
        tempAlpha[count+1] = alpha[levelCount][count+1];
        tempTau[count] = tau[levelCount][count];
        ++count;
      }
      
      // Initialize new quadLevel
      quadLevel myLevel(tempQuad,tempAlpha,tempTau,startIndex);
      
      // Add new quadLevel to quadrature vector
      quadrature.push_back(myLevel);
      
      // Advance iterates
      ++levelCount;	
      startIndex = i;
      myXi = quadSet[i][0];
    
    // The final index in quadSet is a boundary case, which is handled
    // as below
    } else if (i==quadSet.size()-1){
      
      stopIndex = i;
      tempQuad.resize(stopIndex-startIndex+1,vector<double>(4,0.0));
      tempAlpha.resize(stopIndex-startIndex+2,0.0);
      tempTau.resize(stopIndex-startIndex+1,0);
      count = 0;
      tempAlpha[0] = alpha[levelCount][0];
      
      for (int iLevel = startIndex; iLevel < stopIndex + 1; ++iLevel){
      
        tempQuad[count] = quadSet[iLevel];
        tempAlpha[count+1] = alpha[levelCount][count+1];
        tempTau[count] = tau[levelCount][count];
        ++count;
      }
      
      quadLevel myLevel(tempQuad,tempAlpha,tempTau,startIndex);
      quadrature.push_back(myLevel);

    }
  }
}
//==============================================================================

//==============================================================================
/// Calculate number of discrete angles and total weight in quadrature set

void Mesh::calcNumAnglesTotalWeight(){        
  int weightIdx = 3;

  nAngles = 0;
  totalWeight = 0.0;

  for (int iLevel = 0; iLevel < quadrature.size(); ++iLevel){
    nAngles = nAngles + quadrature[iLevel].nOrd;

    for (int iOrd = 0; iOrd < quadrature[iLevel].nOrd; ++iOrd)
      totalWeight = totalWeight + quadrature[iLevel].quad[iOrd][weightIdx];
  }
}
//==============================================================================

//==============================================================================
/// Calculates the ordinate index based on the p and q indices provided
///
/// @param [in] p The first quadrature index
/// @param [in] q The second quadrature index
/// @param [out] index The ordinate index 
int Mesh::quad_index(int p, int q){

  // Index to be returned
  int index; 

  // If the p < 6, xi is negative and we simply call lower quad index
  if(p < 6)
    index = low_quad_index(p,q);

  // Otherwise, xi is positive and we need to perform some algebraic 
  // manipulations to ge the correct index
  else
    index = 42 + 42 - low_quad_index(11 - p, 2*(12 - p) - q);
                  
  return index;
}
//==============================================================================

//=============================================================================
/// Returns an index for an ordinate located where xi < 0 
///
/// @param [in] p The first quadrature index
/// @param [in] q The second quadrature index
/// @param [out] index The ordinate index 
int Mesh::low_quad_index(int p, int q){

  int index = q + p*p + p;
  return index;

}
//==============================================================================

//==============================================================================
/// Prints out quadrature set, differencing coefficients, and tau coefficients

void Mesh::printQuadSet(){        

  // Print quadrature set	
  cout << "QUADRATURE SET:"<<endl;
  cout << "Number of ordinates: " << nAngles;
  cout << ", total weight: " << totalWeight << endl; 
  cout << setw(10) << "xi" << setw(10) <<"mu" << setw(10) << "eta"; 
  cout << setw(10) <<"weight" << setw(10) << "ordIdx" <<endl;

  for (int i = 0; i < quadrature.size(); ++i){
    for(int j = 0; j < quadrature[i].nOrd; ++j){
      for(int k = 0; k < quadrature[i].quad[j].size(); ++k){

        cout << setw(10) <<quadrature[i].quad[j][k];
      } 
    cout<<setw(10)<<quadrature[i].ordIdx[j]<< endl;
    }
  }

  cout << "Number of levels: " << quadrature.size() << endl;
  cout<<""<< endl;

  // Print differencing materials
  cout << "DIFFERENCING COEFFICIENTS:"<<endl;
  for (int i = 0; i < quadrature.size(); ++i){
    for(int j = 0; j < quadrature[i].alpha.size(); ++j){

      cout << setw(12) <<quadrature[i].alpha[j]; 

    }
    cout<< endl;
  }

  // Print tau info	
  cout<< endl;
  cout << "TAU:"<<endl;
  for (int i = 0; i < quadrature.size(); ++i){
    for(int j = 0; j < quadrature[i].tau.size(); ++j){

      cout << setw(12) <<quadrature[i].tau[j]; 

    }
    cout<< endl;
  }
  cout<< endl;

}
//==============================================================================	