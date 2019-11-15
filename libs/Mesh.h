#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <armadillo>
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std;
using namespace arma;

class quadLevel
{
	public:
	quadLevel(vector< vector<double> > myQuad,\
  		vector<double> myAlpha,\
  		vector<double> myTau,\
                int myStartIndex);  	
	//number of ordinates on this quadrature level
  	int nOrd;		
	//quadrature set
  	vector< vector<double> > quad;
        //difference coefficients
  	vector<double> alpha;
	//factor for linear interpolation of half angles
        vector<double> tau;
	//index of each ordinate in angular flux matrices
	vector<int> ordIdx; 
	
};

class Mesh
{
	public:
	Mesh(YAML::Node * myInput);  	
  	int n,nAngles;		
	double dz,dr,Z,R,dt,totalWeight; 
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
 	void calcNumAnglesTotalWeight();
	int quad_index(int p,int q);
	int low_quad_index(int p,int q);
	YAML::Node * input;
};

typedef Eigen::Array<bool,Eigen::Dynamic,1> VectorXb;

#endif
