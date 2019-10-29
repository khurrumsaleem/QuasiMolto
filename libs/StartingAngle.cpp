// File: StartingAngle.cpp
// Purpose: calculate initial half angle fluxes needed for 
// approximation of the angular redistribution term of the
// neutron transport equation in RZ geometry
// Date: October 28, 2019

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "Mesh.h"

using namespace std; 
using namespace arma;

//==============================================================================
//! StartingAngle class that solves RZ neutron transport at the starting angles

class StartingAngle
{
        public:
        // public functions
        StartingAngle(Mesh * myMesh,YAML::Node * myInput);
        void calcStartingAngle();

        private:
        // private functions
        YAML::Node * input;
        Mesh * mesh;
};

//==============================================================================

//==============================================================================
//! StartingAngle object constructor

StartingAngle::StartingAngle(Mesh * myMesh,\
                             YAML::Node * myInput)
{
	// Point to variables for mesh and input file
	mesh = myMesh;
	input = myInput;
};

//==============================================================================

//==============================================================================
//! calcStartingAngle Calculate the starting half angle angular flux

// Solves using simple corner balance on a simplified neutron transport equation 
// in RZ geometry. The simplification comes about by substituting quadrature 
// values that correspond to the starting half angle, and results in the 
// elimination of the angular redistribution term. 

void StartingAngle::calcStartingAngle()
{
        // index xi value is stored in in quadLevel
	const int xiIndex = 0;
        // temporary variable used for looping though quad set
	double xi;
        int zStart;
        int zEnd; 
        int zInc; 
	int borderCellZ;
        int borderCellR;
        vector<int> withinUpstreamR(2);
        vector<int> outUpstreamR(2);
        vector<int> withinUpstreamZ(2);
        vector<int> outUpstreamZ(2);
        double gamma;
	// within cell leakage matrices in R and Z directions
	double kRCoeff;
        double kZCoeff;
	mat kR(4,4,fill::zeros);
        mat kZ(4,4,fill::zeros);
	// out of cell leakage matrices in Rand Z directions
        double lRCoeff;
        double lZCoeff;
	mat lR(4,4,fill::zeros);
        mat lZ(4,4,fill::zeros);
        // reaction matrices
	double t1Coeff;
	double t2Coeff;
	mat t1(4,4,fill::zeros);
        mat t2(4,4,fill::zeros);
	
	for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi){
		// get xi for this quadrature level
		xi = mesh->quadrature[iXi].quad[0][xiIndex];
                rStart = mesh->drs.size()-1;
                borderCellR = 1;
                withinUpstreamR = {1,4};
                outUpstreamR = {2,3};
		// depending on xi, define loop constants	
		if (xi > 0) {
			zStart = mesh->dzs.size()-1;
                        zEnd = 1;
                        zInc = -1;
                        borderCellZ = 1;
                        withinUpstreamZ = {3,4};
                        outUpstreamZ = {1,2};
		}
                else {			
			zStart = 1;
                        zEnd = mesh->dzs.size();
                        zInc = 1;
                        borderCellZ = -1;
                        withinUpstreamZ = {1,2};
                        outUpstreamZ = {3,4};
		}
                for (int iR = mesh->drs.size(); iR > 0; --iR){
                        cout << "iR: " << iR << endl;
			for (int iZ = zStart, count = 0; \
			    count < mesh->dzs.size(); iZ = iZ + zInc, ++count){
                        	cout << "iZ: " << iZ << endl;
				gamma = mesh->rEdge(iR)/mesh->rEdge(iR+1);
				// calculate radial within cell leakage matrix
                                kRCoeff = mesh->dzs(iZ)*mesh->rEdge(iR+1)/8.0;
				calckR();
				// calculate axial within cell leakage matrix
                                kZCoeff = mesh->drs(iR)*mesh->rEdge(iR+1)/16.0;
				calckZ();
				// calculate radial out of cell leakage matrix
                                lRCoeff = mesh->dzs(iZ)*mesh->rEdge(ir+1)/2.0;
				calclR();
				// calculate axial out of cell leakage matrix
                                lZCoeff = mesh->drs(iR)*mesh->rEdge(ir+1)/8.0;
				calclZ();
				// calculate first collision matrix
                                t1Coeff = mesh->drs(iR)*mesh->dzs(iZ)*mesh->rEdge(ir+1)/16.0;
				calct1();
				// calculate second collision matrix
                                t2Coeff = mesh->drs(iR)*mesh->dzs(iZ)/4.0;
				calct2();
                        }
		}
	}
	

};

//==============================================================================
