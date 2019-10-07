
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class Mesh
{
	public:
	//order of quadrature set
  	int n;		
	//quadrature set
  	vector<double> quadSet;
	//mu variable in quadrature et
	vector<double> mu;
	//ordi
  	vector< vector<double> > ordinates;
	Mesh();  	
  	void calcMu();	
  	void calcQuadSet();
	void printQuadSet();
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
        
void Mesh::calcQuadSet(){
	calcMu();
	// Find number of sets of ordinates magnitude is one
	int counter = 0;
	for (int i = 0; i < n/2; ++i){
		for (int j =0; j < n/2; ++j){
			for (int k =0; k < n/2; ++k){
				if (abs(pow(mu[i],2.0)+pow(mu[j],2.0)+pow(mu[k],2.0)-1)<1E-5){
				++counter;
				}
			}
		}
	}

	ordinates.resize(counter,vector<double>(4,0.0)); 
	counter = 0;
	vector<double> thisRow;
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

	double weight1 = 0.0707626;
	double weight2 = 0.0558811;
	double weight3 = 0.0373377;
	double weight4 = 0.0502819;
	double weight5 = 0.0258513;
	vector<int> ord1 = {1,6,21}; 
	vector<int> ord2 = {2,5,7,11,19,20};
	vector<int> ord3 = {3,4,12,15,16,18};
	vector<int> ord4 = {8,10,17};
	vector<int> ord5 = {9,13,14};
	for(int i = 0; i < ord1.size(); ++i){
		ordinates[ord1[i]-1][3] = weight1;
	 }
	 for(int i = 0; i < ord2.size(); ++i)
		 ordinates[ord2[i]][3] = weight2;
 
	 for(int i = 0; i < ord3.size(); ++i)
		 ordinates[ord3[i]][3] = weight3;
 
	 for(int i = 0; i < ord4.size(); ++i)
		 ordinates[ord4[i]][3] = weight4;
 
	 for(int i = 0; i < ord5.size(); ++i)
		 ordinates[ord5[i]][3] = weight5;

}

void Mesh::printQuadSet(){        	
	// Check to make sure mu is initialized before printing
	cout << "=====QUADRATE SET EDITS=====" << endl;
	if (mu.size()>0){
		for(int imu = 0; imu < n/2; ++imu){
        		cout << "mu(" << imu << ")=" << mu[imu] << endl;
		}
	} else {
        cout << "ERROR: mu of quadrature set not initialized! " << endl;
	}
	cout << ordinates.size() << endl;
}	
