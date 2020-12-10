#ifndef MESH_H
#define MESH_H
#define EIGEN_SUPERLU_SUPPORT

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>
#include <armadillo>
#include <petsc.h>
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"
#include "../TPLs/eigen-git-mirror/Eigen/IterativeLinearSolvers"
#include "../TPLs/eigen-git-mirror/unsupported/Eigen/IterativeSolvers"
#include "../TPLs/eigen-git-mirror/Eigen/SuperLUSupport"

using namespace std;

class WriteData; // forward declaration

class qdCell
{

	public:
	qdCell();  	
	//number of ordinates on this quadrature level
  	int cFIndex=0,nFIndex=0,sFIndex=0,eFIndex=0,wFIndex=0;
        int nCIndex=0,sCIndex=0,eCIndex=0,wCIndex=0;

};


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
  	int n,nAngles,nR,nZ;		
        int state = 1; 
        double dz,dr,drCorner,dzCorner,Z,R,dt,T,totalWeight; 
        bool verbose = false,petsc = false;
  	vector< vector<double> > quadSet;
  	vector< vector<double> > alpha;
        vector< vector<double> > tau;
        vector< vector<double> > cellVol;
        vector< vector<double> > cellVSA;
        vector<double> dts;
        vector<double> ts;
        vector<bool> outputOnStep;
        arma::rowvec dzs,dzsCorner;
        arma::rowvec drs,drsCorner;
	arma::rowvec rEdge;
	arma::rowvec zEdge;
        arma::rowvec rCent;
        arma::rowvec zCent;
	arma::rowvec rCornerEdge;
	arma::rowvec zCornerEdge;
        arma::rowvec rCornerCent;
        arma::rowvec zCornerCent;
        arma::rowvec rVWCornerCent;
        vector<quadLevel> quadrature;
	vector<qdCell> qdCells; 
        Eigen::MatrixXd volume;
        string outputDir = "mesh/";
        WriteData * output;

        // Recirculation loop parameters
        double dzCornerRecirc,recircZ;
        int nZrecirc;
	arma::rowvec zCornerEdgeRecirc,dzsCornerRecirc;
        arma::rowvec zCornerCentRecirc;

        // Functions
  	vector<int> getQDCellIndices(int iR, int iZ);
  	vector<double> getGeoParams(int iR, int iZ);
  	vector<double> getRecircGeoParams(int iR, int iZ);
  	void advanceOneTimeStep();
  	void calcQuadSet();
        void writeVars();
	void printQuadSet();
        void checkOptionalParams();
	
        private:
	vector<double> mu;
  	vector< vector<double> > ordinates;
  	void calcMu();	
	void calcAlpha();
        void calcTau();
        void calcSpatialMesh();
        void calcQDCellIndices(int nCornersR,int nCornersZ);
        void addLevels();
 	void calcNumAnglesTotalWeight();
 	void calcTimeMesh();
  	void calcRecircMesh();
	int quad_index(int p,int q);
	int low_quad_index(int p,int q);
	int getQDCellIndex(int iR,int iZ);
	YAML::Node * input;
};

typedef Eigen::Array<bool,Eigen::Dynamic,1> VectorXb;

#endif
