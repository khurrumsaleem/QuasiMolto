#ifndef MATERIAL_H
#define MATERIAL_H

#include "Mesh.h"

using namespace std; 

//==============================================================================
//! Material class that holds material and geometry information

class Material
{
        public:
	int matID;
	string name;
 	Eigen::VectorXd chiP,chiD;
        double density; // density of the material
        double gamma; // fraction of energy deposited from gamma rays
        double k; // thermal conductivity
        double cP; // specific heat
        double omega; // energy release per fission	
        bool stationary = true; // is the material stationary
        
        // public functions
	Material(int myMatID,\
		string name,\
                vector<Eigen::MatrixXd> mySigT,\
                vector<Eigen::MatrixXd> mySigS,\
                vector<Eigen::MatrixXd> mySigF,\
                vector<Eigen::MatrixXd> myNu,\
                vector<Eigen::MatrixXd> myNeutV,\
                Eigen::MatrixXd myFlowVelocity,\
                Eigen::VectorXd myChiP,\
                Eigen::VectorXd myChiD,\
                double myDensity,\
                double myGamma,\
                double myK,\
                double mycP,\
                double myOmega,\
                bool myStationary=true);
        double getSigT(int eIdx,double temp);
        double getSigF(int eIdx,double temp);
        double getSigS(int eIdxPrime,int eIdx,double temp);
        double getNu(int eIdx,double temp);
        double getNeutV(int eIdx,double temp);
        double getFlowVelocity(double time);
        double interpolateParameter(Eigen::MatrixXd param,double temp);
        void checkMat();
	void edit();

        private:
        vector<Eigen::MatrixXd> sigT,sigF,sigS,neutV,nu;
        Eigen::MatrixXd flowVelocity;

};

//==============================================================================

#endif
