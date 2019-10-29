#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <armadillo>
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"

using namespace std;
using namespace arma;

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
	Mesh(YAML::Node * myInput);  	
  	int n;		
	double dz; 
	double dr;
	double Z; 
	double R;
  	vector< vector<double> > quadSet;
  	vector< vector<double> > alpha;
        vector< vector<double> > tau;
        vector< vector<double> > cellVol;
        vector< vector<double> > cellVSA;
        rowvec dzs;
        rowvec drs;
	rowvec rEdge;
	rowvec zEdge;
        rowvec rCent;
        rowvec zCent;
        vector<quadLevel> quadrature;
  	void calcQuadSet();
	void printQuadSet();
	
        private:
	vector<double> mu;
  	vector< vector<double> > ordinates;
  	void calcMu();	
	void calcAlpha();
        void calcTau();
        void calcSpatialMesh();
        void addLevels();
	int quad_index(int p,int q);
	int low_quad_index(int p,int q);
	YAML::Node * input;
};


#endif
