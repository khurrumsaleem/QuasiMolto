
#include <iostream>

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
		double constant = 2*(1-3*mu[1]^2)/(n-2);
		for (int imu = 1; i < n; ++i){
			mu[imu] = sqrt(mu[imu - 1]^2 + c);
		}
		cout << "mu(:) :"<<mu<<endl;
}
        
void Mesh::calcQuadSet(){
	calcMu();
}

void Mesh::printQuadSet(){        	
	cout << "mu(:) :"<<mu<<endl;
}	
