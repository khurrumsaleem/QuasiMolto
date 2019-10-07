
#include <iostream>
#include <cmath>

using namespace std;

class Mesh
{
	public:
	//order of quadrature set
  	int n;		
	//quadrature set
  	double* quadSet;
	//mu variable in quadrature et
  	double* mu;
	//ordi
  	double* ordinates;
	Mesh();  	
  	void calcMu();	
  	void calcQuadSet();
	void printQuadSet();
};
// Contructor for this object. 
Mesh::Mesh(){
	n=12; // currently only works for a quadrature set of order 12
	quadSet = 0; //initialize to null
  	mu = 0; //initialize to null
  	ordinates = 0; //initialize to null
}

void Mesh::calcMu(){
		// Allocate memory and fix one degree of freedom
    		mu = new double[n/2];
		mu[0] = 0.1672126;
		double myConstant = 2.0*(1.0-3.0*pow(mu[0],2.0))/(n-2.0);
		for (int imu = 1; imu < n/2; ++imu){
			mu[imu] = sqrt(pow(mu[imu - 1],2.0) + myConstant);
		}
}
        
void Mesh::calcQuadSet(){
	calcMu();
}

void Mesh::printQuadSet(){        	
	// Check to make sure mu is initialized before printing
	cout << "=====QUADRATE SET EDITS=====" << endl;
	if(mu){
		for(int imu = 0; imu < n/2; ++imu){
        		cout << "mu(" << imu << ")=" << mu[imu] << endl;
		}
	} else {
        cout << "ERROR: mu of quadrature set not initialized! " << endl;
	}
}	
