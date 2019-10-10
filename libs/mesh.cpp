// File: mesh.cpp
// Purpose: defines a class to hold mesh data such as spatial mesh and 
// quadrature set
// Author: Aaron James Reynolds
// Date: October 9, 2019

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;
//==============================================================================
//! Mesh class that contains mesh data for a simulation

class Mesh
{
	public:
	Mesh();  	
	//order of quadrature set
  	int n;		
	//quadrature set
  	vector< vector<double> > quadSet;
        //difference coefficients
  	vector< vector<double> > alpha;
        // public functions
  	void calcQuadSet();
	void printQuadSet();
	
        private:
	// viable ordinates
	vector<double> mu;
	// all three ordinates and weights
  	vector< vector<double> > ordinates;
        // private functions
  	void calcMu();	
	void calcAlpha();
	int quad_index(int p,int q);
	int low_quad_index(int p,int q);
};
//==============================================================================

//==============================================================================
//! Mesh Contructor for Mesh object. 
Mesh::Mesh(){
	n=12; // currently only works for a quadrature set of order 12 
}
//==============================================================================

//==============================================================================
//! calcQuadSet function that calculates a quadrature set and 
//! differencing coefficients for use in RZ geometry

//! Uses the methodology laid out Methods of Computational Transport by
//! Lewis and Miller. A level symmetric quadrature set is calculated. The p,q 
//! indexing scheme is defined on page 166 in Figure 4-7. Differencing 
//! coefficients are then calculated using that quadrature set.
//! TODO: a lot of this is very specific to an n=12 quadrature set. Would be 
//! better if it was more general          
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
	int counter;
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
	tempQuad.resize(4*counter,vector<double>(3,0.0)); 
        for (int iquad = 0; iquad < 4; ++iquad){
		for (int iset = 0; iset < ordinates.size(); ++iset){
			tempRow = {coeffQ1[iquad]*ordinates[iset][0],\
				coeffQ2[iquad]*ordinates[iset][1],\
			        ordinates[iset][3]};
			tempQuad[counter*iquad+iset] = tempRow;
		}
	}

        // Sort the quadrature set in ascending order	
	orderedQuad.resize(4*counter,vector<double>(3,0.0)); 
	int rowIndexOfMinimum;
	for (int iOrdered = 0; iOrdered < orderedQuad.size(); ++iOrdered){
		rowIndexOfMinimum = 0;
		for (int iTemp = 0; iTemp < tempQuad.size(); ++iTemp){
			if (tempQuad[iTemp][0]<tempQuad[rowIndexOfMinimum][0]){
				rowIndexOfMinimum = iTemp;
			} else if (tempQuad[iTemp][0]==tempQuad[rowIndexOfMinimum][0]){
				if (tempQuad[iTemp][1]<tempQuad[rowIndexOfMinimum][1]){
					rowIndexOfMinimum = iTemp;	
				}
			}
		}
        	// Place the next smallest on orderedQuad
		orderedQuad[iOrdered] = tempQuad[rowIndexOfMinimum];
        	// Set this last placed row to some large values so we don't
        	// accidentally recognize it as a minimum again
		tempQuad[rowIndexOfMinimum] = {100,100,100};
	}
        // Set quadSet
	quadSet.resize(4*counter,vector<double>(3,0.0)); 
	quadSet = orderedQuad;
        // Calculate differencing coefficients knowing the quadSet
	calcAlpha();
        
}
//==============================================================================

//==============================================================================
//! calcMu calculate suitable ordinates

//! 
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
//! calcAlpha function for calculating differencing coefficients in RZ geometry

//! based on approach in Lewis and Miller
void Mesh::calcAlpha(){
	
	vector<int> rowLength = {2,4,6,8,10,12,12,10,8,6,4,2};
	alpha.resize(12,vector<double>(13,0.0)); 
	
        for (int i = 0; i < alpha.size(); ++i){
		alpha[i][0] = 0; // initialize to 0 on the edge case to conserve
                                 // neutrons. Explanation on page 179 of L&M.
		for (int j = 0; j < rowLength[i]; ++j){
			alpha[i][j+1] = alpha[i][j]-quadSet[quad_index(i,j)][1]\
					*quadSet[quad_index(i,j)][2];
			// force sufficiently small quantities to 0
                        // keeps negative values resulting from finite 
                        // precision from coming up, too
                        if (alpha[i][j+1]<1E-15)
                          alpha[i][j+1] = 0.0;
		}
	}
}
//==============================================================================

//==============================================================================
//! quad_index returns a sequential index number based on the p and q indices 
//! provided as arguments

//! \param p the first quadrature index
//! \param q the second quadrature index
int Mesh::quad_index(int p, int q){
	int index; ///< index to be returned
        // if the p < 6, xi is negative and we simply call lower quad index
	if(p < 6 )
		index = low_quad_index(p,q);
        // otherwise, xi is positive and we need to perform some algebraic 
        // manipulations to ge the correct index
	else
		index = 42 + 42 - low_quad_index(11 - p, 2*(12 - p) - q);
			
	return index;
}
//==============================================================================

//=============================================================================
//! low_quad_index returns an index considering sequential numbering in the 
//! region where xi is negative

//! \param p the first quadrature index
//! \param q the second quadrature index
int Mesh::low_quad_index(int p, int q){
	int index = q + p*p + p;
	return index;
}
//==============================================================================

//==============================================================================
//! editAngularMesh prints out quadrature set and differencing coefficients
void Mesh::printQuadSet(){        	
	cout << "QUADRATURE SET:"<<endl;
	cout << setw(10) << "mu" << setw(10) <<"xi" << setw(10) <<"weight" << endl;
	for (int i = 0; i < quadSet.size(); ++i){
		for(int j = 0; j < quadSet[i].size(); ++j){
        		cout << setw(10) <<quadSet[i][j]; 
		}
		cout<<""<< endl;
	}
	cout << "Size: " << quadSet.size() << "x" << quadSet[0].size() << endl;
	cout<<""<< endl;
	cout << "DIFFERENCING COEFFICIENTS:"<<endl;
	for (int i = 0; i < alpha.size(); ++i){
		for(int j = 0; j < alpha[i].size(); ++j){
        		cout << setw(12) <<alpha[i][j]; 
		}
		cout<<""<< endl;
	}
	cout << "Size: " << alpha.size() << "x" << alpha[0].size() << endl;

}
//==============================================================================	
