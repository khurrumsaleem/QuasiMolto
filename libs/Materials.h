#ifndef MATERIALS_H
#define MATERIALS_H

#include "Mesh.h"
#include "Material.h"
#include "CollapsedCrossSections.h"

using namespace std; 

//==============================================================================
//! Material class that holds material and geometry information

class Materials
{
        public:
        CollapsedCrossSections * oneGroupXS;
	Eigen::MatrixXi matMap; //resized later
	Eigen::MatrixXd flowVelocity; 
	Eigen::MatrixXd temperature; 
	Eigen::MatrixXd recircFlowVelocity; 
        Eigen::VectorXd neutV;
        int nGroups;
        bool posVelocity;
        bool uniformTemperature = false;
        double uniformTempValue;
        // public functions
        Materials(Mesh * myMesh,YAML::Node * myInput);
	void readMats();
	void readGeom();
	void readFlowVelocity(double time = 0);
        void setMatRegion(int myIndex,double rIn,double rOut,\
		double zUp,double zLow);
	double sigT(int zIdx,int rIdx,int eIndx);
	double zSigT(int zIdx,int rIdx,int eIndx);
	double rSigT(int zIdx,int rIdx,int eIndx);
	double sigS(int zIdx,int rIdx,int gprime,int g);
	double sigF(int zIdx,int rIdx,int eIndx);
	double chiP(int zIdx,int rIdx,int eIndx);
	double chiD(int zIdx,int rIdx,int eIndx);
	double nu(int zIdx,int rIdx,int eIdx);
	double neutVel(int zIdx,int rIdx,int eIdx);
	double zNeutVel(int zIdx,int rIdx,int eIdx);
	double rNeutVel(int zIdx,int rIdx,int eIdx);
	double density(int zIdx,int rIdx);
	double gamma(int zIdx,int rIdx);
	double k(int zIdx,int rIdx);
	double cP(int zIdx,int rIdx);
	double omega(int zIdx,int rIdx);
	double coreFlowVelocity(int zIdx,int rIdx, double time = 0);
        void updateTemperature(Eigen::MatrixXd myTemp);
        vector<Eigen::MatrixXd> readTempDependentYaml(string file);
        Eigen::MatrixXd readTimeDependentYaml(string file);
        void checkMats();
        void initCollapsedXS();
        void edit();
	vector<shared_ptr<Material>> matBank;

        private:
        // private variables
        YAML::Node * input;
        Mesh * mesh;
	map<string,int> mat2idx;

};

//==============================================================================

#endif
