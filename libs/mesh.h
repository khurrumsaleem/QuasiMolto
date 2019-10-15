#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>

using namespace std;

class quadLevel
{
	public:
	quadLevel(vector< vector<double> > myQuad,\
  		vector<double> myAlpha,\
  		vector<double> myTau);  	
	//number of ordinates on this quadrature level
  	int nOrd;		
	//quadrature set
  	vector< vector<double> > quad;
        //difference coefficients
  	vector<double> alpha;
	//factor for linear interpolation of half angles
        vector<double> tau;
	
};

class Mesh
{
	public:
	Mesh();  	
  	int n;		
  	vector< vector<double> > quadSet;
  	vector< vector<double> > alpha;
        vector< vector<double> > tau;
        vector<quadLevel> quadrature;
  	void calcQuadSet();
	void printQuadSet();
	
        private:
	vector<double> mu;
  	vector< vector<double> > ordinates;
  	void calcMu();	
	void calcAlpha();
        void calcTau();
        void addLevels();
	int quad_index(int p,int q);
	int low_quad_index(int p,int q);
};


#endif
