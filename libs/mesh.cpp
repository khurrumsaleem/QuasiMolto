
#include <iostream>
#include <cmath>

using namespace std;

class Mesh
{
	public:
	
  	int n=12;		
  	double* quadSet = NULL;
  	double* mu = NULL;
  	double* ordinates = NULL;
  	
  	void calcMu();	
  	void calcQuadSet();
	void printQuadSet();
};

void Mesh::calcMu(){
		// Allocate memory and fix one degree of freedom
    		mu = new double[n/2];
		mu[0] = 0.1672126;
		double myConstant = 2.0*pow(1.0-3.0*mu[1],2.0)/(n-2.0);
		for (int imu = 1; imu < n/2; ++imu){
			mu[imu] = sqrt(pow(mu[imu - 1],2.0) + myConstant);
			cout << "mu("<<imu<<"):"<<mu[imu]<<endl;
		}
		cout << "mu(:) :"<<mu<<endl;
}
        
void Mesh::calcQuadSet(){
	calcMu();
}

void Mesh::printQuadSet(){        	
	cout << "mu(:) :"<<mu<<endl;
}	
