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
        mat calckR(double myGamma);
        mat calckZ(double myGamma);
        mat calclR(double myGamma);
        mat calclZ(double myGamma);
        mat calct1(double myGamma);
        mat calct2(double myGamma);
        rowvec calcSubCellVol(int myiZ, int myiR);

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
	double sqrtXi;
        int zStart;
        int rStart;
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
	// A matrix of linear system Ax=b
	mat A(4,4,fill::zeros);
	rowvec subCellVol;

	// need to bring this in from materials... kluge for now
        double sigT = 10.0;
	
	for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi){
		// get xi for this quadrature level
		xi = mesh->quadrature[iXi].quad[0][xiIndex];
                sqrtXi = pow(1-pow(xi,2),0.5); 
                rStart = mesh->drs.size()-1;
                borderCellR = 1;
                withinUpstreamR = {1,4};
                outUpstreamR = {2,3};
		// depending on xi, define loop constants	
		if (xi > 0) {
			zStart = mesh->dzs.size()-1;
                        zEnd = 0;
                        zInc = -1;
                        borderCellZ = 1;
                        withinUpstreamZ = {3,4};
                        outUpstreamZ = {1,2};
		}
                else {			
			zStart = 0;
                        zEnd = mesh->dzs.size();
                        zInc = 1;
                        borderCellZ = -1;
                        withinUpstreamZ = {1,2};
                        outUpstreamZ = {3,4};
		}
                for (int iR = rStart; iR > 0; --iR){
                        cout << "iR: " << iR << endl;
			for (int iZ = zStart, count = 0; \
			    count < mesh->dzs.size(); iZ = iZ + zInc, ++count){
                        	cout << "iZ: " << iZ << endl;
				gamma = mesh->rEdge(iR)/mesh->rEdge(iR+1);
				// calculate radial within cell leakage matrix
                                kRCoeff = mesh->dzs(iZ)*mesh->rEdge(iR+1)/8.0;
				kR=calckR(gamma);
				kR=kRCoeff*kR;
				// calculate axial within cell leakage matrix
                                kZCoeff = mesh->drs(iR)*mesh->rEdge(iR+1)/16.0;
				kZ=calckZ(gamma);
				kZ=kZCoeff*kZ;
				// calculate radial out of cell leakage matrix
                                lRCoeff = mesh->dzs(iZ)*mesh->rEdge(iR+1)/2.0;
				lR=calclR(gamma);
				lR=lRCoeff*lR;
				// calculate axial out of cell leakage matrix
                                lZCoeff = mesh->drs(iR)*mesh->rEdge(iR+1)/8.0;
				lZ=calclZ(gamma);
				lZ=lZCoeff*lZ;
				// calculate first collision matrix
                                t1Coeff = mesh->drs(iR)*mesh->dzs(iZ)*mesh->rEdge(iR+1)/16.0;
				t1=calct1(gamma);
				t1=t1Coeff*t1;
				// calculate second collision matrix
                                t2Coeff = mesh->drs(iR)*mesh->dzs(iZ)/4.0;
				t2=calct2(gamma);
				t2=t2Coeff*t2;
              			
				A = sqrtXi*kR+xi*kZ+sigT*t1+sqrtXi*t2;
				subCellVol = calcSubCellVol(iZ,iR);	
                        }
		}
	}
	

};

mat StartingAngle::calckR(double myGamma){
	double a = -(1+myGamma);
	double b = 1+myGamma;
	mat kR(4,4,fill::zeros);

	kR(0,0) = a; kR(0,1) = a;
	kR(1,0) = b; kR(1,1) = b;
	kR(2,2) = b; kR(2,3) = b;
	kR(3,2) = a; kR(3,3) = a;
         
	return kR;
}

mat StartingAngle::calckZ(double myGamma){
	double a = 1+3*myGamma;
	double b = 3+myGamma;
	mat kZ(4,4,fill::zeros);

	kZ(0,0) = a; kZ(0,3) = a;
	kZ(1,2) = b; kZ(1,3) = b;
	kZ(2,2) = -b; kZ(2,3) = -b;
	kZ(3,0) = -a; kZ(3,3) = -a;
         
	return kZ;
}

mat StartingAngle::calclR(double myGamma){
	double a = myGamma;
	double b = -1;
	mat lR(4,4,fill::zeros);

	lR(0,0) = a; 
	lR(1,1) = b; 
	lR(2,2) = b; 
	lR(3,3) = a; 
         
	return lR;
}

mat StartingAngle::calclZ(double myGamma){
	double a = 1+3*myGamma;
	double b = 3+myGamma;
	mat lZ(4,4,fill::zeros);

	lZ(0,0) = -a; 
	lZ(1,1) = -b; 
	lZ(2,2) = b; 
	lZ(3,3) = a; 
         
	return lZ;
}

mat StartingAngle::calct1(double myGamma){
	double a = 1+3*myGamma;
	double b = 3+myGamma;
	mat t1(4,4,fill::zeros);

	t1(0,0) = a; 
	t1(1,1) = b; 
	t1(2,2) = b; 
	t1(3,3) = a; 
         
	return t1;
}

mat StartingAngle::calct2(double myGamma){
	double a = 1;
	mat t2(4,4,fill::zeros);

	t2(0,0) = a; 
	t2(1,1) = a; 
	t2(2,2) = a;
	t2(3,3) = a; 
         
	return t2;
}

rowvec StartingAngle::calcSubCellVol(int myiZ, int myiR){
	rowvec subCellVol(3,fill::zeros);
	
	subCellVol(0) = (mesh->dzs(myiZ)/2)*(pow(mesh->rCent(myiR),2)-\
		pow(mesh->rEdge(myiR),2));
	subCellVol(1) = (mesh->dzs(myiZ)/2)*(pow(mesh->rEdge(myiR+1),2)-\
                pow(mesh->rCent(myiR),2));
	subCellVol(2) = (mesh->dzs(myiZ)/2)*(pow(mesh->rEdge(myiR+1),2)-\
                pow(mesh->rCent(myiR),2));
	subCellVol(0) = (mesh->dzs(myiZ)/2)*(pow(mesh->rCent(myiR),2)-\
		pow(mesh->rEdge(myiR),2));

	return subCellVol;
}
//==============================================================================
