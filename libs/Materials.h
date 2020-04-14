#ifndef MATERIALS_H
#define MATERIALS_H

#include "Mesh.h"
#include "Material.h"

using namespace std; 
using namespace arma;

//==============================================================================
//! Material class that holds material and geometry information

class Materials
{
        public:
	Eigen::MatrixXi matMap; //resized later
	Eigen::MatrixXd flowVelocity; 
        Eigen::VectorXd neutV;
        int nGroups;
        bool posVelocity;
        // public functions
        Materials(Mesh * myMesh,YAML::Node * myInput);
	void readMats();
	void readGeom();
	void readFlowVelocity();
        void setMatRegion(int myIndex,double rIn,double rOut,\
		double zUp,double zLow);
	double sigT(int zIdx,int rIdx,int eIndx);
	double sigS(int zIdx,int rIdx,int gprime,int g);
	double sigF(int zIdx,int rIdx,int eIndx);
	double chiP(int zIdx,int rIdx,int eIndx);
	double chiD(int zIdx,int rIdx,int eIndx);
	double nu(int zIdx,int rIdx);
	double density(int zIdx,int rIdx);
	double gamma(int zIdx,int rIdx);
	double k(int zIdx,int rIdx);
	double cP(int zIdx,int rIdx);
	double omega(int zIdx,int rIdx);
        void checkMats();
        void edit();

        private:
        // private functions
        YAML::Node * input;
        Mesh * mesh;
	map<string,int> mat2idx;
	vector<shared_ptr<Material>> matBank;

};

//==============================================================================

#endif
