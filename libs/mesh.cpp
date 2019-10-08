
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

class Mesh
{
	public:
	//order of quadrature set
  	int n;		
	//quadrature set
  	vector< vector<double> > quadSet;
	//mu variable in quadrature et
	vector<double> mu;
	//ordi
  	vector< vector<double> > ordinates;
  	vector< vector<double> > alpha;
	Mesh();  	
  	void calcMu();	
  	void calcQuadSet();
	void calcAlpha();
	void printQuadSet();
	int quad_index(int p,int q);
	int low_quad_index(int p,int q);
};

// Contructor for this object. 
Mesh::Mesh(){
	n=12; // currently only works for a quadrature set of order 
}

void Mesh::calcMu(){
		// Allocate memory and fix one degree of freedom
    		mu.resize(n/2,0.0);
		mu[0] = 0.1672126;
		double myConstant = 2.0*(1.0-3.0*pow(mu[0],2.0))/(n-2.0);
		for (int imu = 1; imu < n/2; ++imu){
			mu[imu] = sqrt(pow(mu[imu - 1],2.0) + myConstant);
		}
}
//
//=============================================================================
//! calcQuadSet Function that calculates a quadrature set and 
//! differencing coefficients for use in RZ geometry

//! Uses the methodology laid out Lewis and Miller. A level symmetric
//! quadrature set is calculated.
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
	orderedQuad[iOrdered] = tempQuad[rowIndexOfMinimum];
	tempQuad[rowIndexOfMinimum] = {100,100,100};
	}
	quadSet.resize(4*counter,vector<double>(3,0.0)); 
	quadSet = orderedQuad;
	calcAlpha();
        
}
//! calcAlpha function for calculating differencing coefficients
void Mesh::calcAlpha(){
	
	vector<int> rowLength = {2,4,6,8,10,12,12,10,8,6,4,2};
	alpha.resize(12,vector<double>(13,0.0)); 
        
        cout << quad_index(7,0) << endl;	
        cout << quad_index(7,1) << endl;	
        cout << quad_index(8,3) << endl;	
        cout << quad_index(0,0) << endl;	
        cout << quad_index(0,1) << endl;	
	for (int i = 0; i < alpha.size(); ++i){
		alpha[i][0] = 0;
		for (int j = 0; j < rowLength[i]; ++j){
			alpha[i][j+1] = alpha[i][j]-quadSet[quad_index(i,j)][1]\
					*quadSet[quad_index(i,j)][2];
		}
	}
}

int Mesh::quad_index(int p, int q){
	int index;
	if(p < 6 )
		index = low_quad_index(p,q);
	else
		index = 42 + 42 - low_quad_index(11 - p, 2*(12 - p) - q);
			
	return index;
}

int Mesh::low_quad_index(int p, int q){
	int index = q + p*p + p;
	return index;
}
//=============================================================================
void Mesh::printQuadSet(){        	
	// Check to make sure mu is initialized before printing
	cout << "QUADRATURE SET:"<<endl;
	cout << setw(10) << "mu" << setw(10) <<"xi" << setw(10) <<"weight" << endl;
	for (int i = 0; i < quadSet.size(); ++i){
		for(int j = 0; j < quadSet[i].size(); ++j){
        		cout << setw(10) <<quadSet[i][j]; 
		}
		cout<<""<< endl;
	}
	
	cout << setw(10) << "mu" << setw(10) <<"xi" << setw(10) <<"weight" << endl;
	for (int i = 0; i < alpha.size(); ++i){
		for(int j = 0; j < alpha[i].size(); ++j){
        		cout << setw(10) <<alpha[i][j]; 
		}
		cout<<""<< endl;
	}
	cout << "Number of entries: " << alpha.size() << endl;

}	
