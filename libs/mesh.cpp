
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

	ordinates.resize(counter,vector<double>(3,0.0)); 
	counter = 0;
	vector<double> thisRow;
	for (int i = 0; i < n/2; ++i){
		for (int j = 0; j < n/2; ++j){
			for (int k = 0; k < n/2; ++k){
				if (abs(pow(mu[i],2.0)+pow(mu[j],2.0)+pow(mu[k],2.0)-1)<1E-5){
				thisRow = {mu[i],mu[j],mu[k]};
				ordinates[counter] = thisRow;
				++counter;
				}
			}
		}
	}

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
