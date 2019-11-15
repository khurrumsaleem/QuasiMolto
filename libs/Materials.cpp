// File: Materials.cpp
// Purpose: contain material characteristics and geometry for a simulation
// Date: October 28, 2019

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "Mesh.h"
#include "Material.h"
#include "Materials.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std; 
using namespace arma;

//==============================================================================
//! Material class object constructor

Materials::Materials(Mesh * myMesh,\
                             YAML::Node * myInput)
{
        vector<double> neutVInp;    

	// Point to variables for mesh and input file
	mesh = myMesh;
	input = myInput;
	matMap.setZero(mesh->zCent.size(),mesh->rCent.size());
	readMats();
	readGeom();
        edit();
        checkMats();        
       
        // check if neutron velocities are specified in input
        if ((*input)["parameters"]["neutron velocity"]){
	        neutVInp = (*input)["parameters"]["neutron velocity"]\
                        .as<vector<double>>();
                neutV.setZero(nGroups);
                for (int iGroup = 0; iGroup < nGroups; ++iGroup){
                        neutV(iGroup) = neutVInp[iGroup];
                }
        } else {
          
                // set default neutron velocities 
                neutV.setOnes(nGroups);
                neutV = 2200.0*neutV;
        }
};

//==============================================================================

//==============================================================================
//! readInput Parse input and define  geometry parameters

void Materials::readGeom()
{
	YAML::Node geometry = (*input)["geometry"];
	string regType,matName;
	double rIn,rOut,zLow,zUp;
	int index;
        
	for (YAML::const_iterator it=geometry.begin(); it!=geometry.end(); ++it){
		regType=it->first.as<string>();
		if (regType=="background"){
			index = mat2idx[it->second.as<string>()];
			setMatRegion(index,0.0,mesh->R,0.0,mesh->Z);
		} else if (regType=="region"){
			index = mat2idx[it->second["material"].as<string>()];
			rIn = it->second["inner-r"].as<double>(); 
			rOut = it->second["outer-r"].as<double>(); 
			zLow = it->second["lower-z"].as<double>(); 
			zUp = it->second["upper-z"].as<double>();
			setMatRegion(index,rIn,rOut,zLow,zUp);
		}
	}
};

//==============================================================================

//==============================================================================
//! readInput Parse input and define material parameters

void Materials::readMats()
{
	YAML::Node mats = (*input)["materials"];
	string name;
	vector<double> sigTInp,sigSInp,sigFInp,chiPInp,chiDInp;
	Eigen::VectorXd sigT,sigF,chiP,chiD; 
	Eigen::MatrixXd sigS;
	int ID,size;
	double nu;
	
        int iCount=0;
	for (YAML::const_iterator it=mats.begin(); it!=mats.end(); ++it){
		mat2idx.insert(pair<string,int>(it->first.as<string>(),iCount));
		name = it->first.as<string>();
		// pull values from input
		sigTInp = it->second["sigT"].as<vector<double>>();
		sigSInp = it->second["sigS"].as<vector<double>>();
		sigFInp = it->second["sigF"].as<vector<double>>();
		chiPInp = it->second["chiP"].as<vector<double>>();
		chiDInp = it->second["chiD"].as<vector<double>>();
		nu = it->second["nu"].as<double>();
		// set size of arma vectors
		size = sigTInp.size();
		sigT.setZero(size); 
		sigF.setZero(size);
		chiP.setZero(size);
		chiD.setZero(size);
		sigS.setZero(size,size);
		// load standard vector inputs into arma vector
		for (int iSig = 0; iSig < size; ++iSig){
			sigT(iSig) = sigTInp[iSig];
			sigF(iSig) = sigFInp[iSig];
                        chiP(iSig) = chiPInp[iSig];
                        chiD(iSig) = chiDInp[iSig];
			for(int iGroup = 0; iGroup < size; ++iGroup){
				sigS(iSig,iGroup) = sigSInp[iSig*size+iGroup];
			}
		}
		// add material to bank
		shared_ptr<Material> newMat (new Material(iCount,name,sigT,sigS,\
                  sigF,chiP,chiD,nu));
		matBank.push_back(std::move(newMat));
		++iCount;
	}	
        nGroups = size;
};

//==============================================================================

//==============================================================================
//! setMatRegions set material region to provided index

// Set all locations within the input bounds equal to the input index
void Materials::setMatRegion(int myIndex,double rIn,double rOut,double zLow,double zUp)
{
	for (int iZ = 0; iZ < mesh->zCent.size(); iZ++){
		if (mesh->zCent(iZ)>zLow && mesh->zCent(iZ)<zUp){
			for (int iR = 0; iR < mesh->rCent.size(); iR++){
				if (mesh->rCent(iR)>rIn && mesh->rCent(iR)<rOut){
					matMap(iZ,iR) = myIndex;
				}
			}
		}
	}
};

//==============================================================================

//==============================================================================
//! sigT return total cross section

// return total cross section at location indicated by input indices
double Materials::sigT(int zIdx,int rIdx,int eIndx){

	double sigT = matBank[matMap(zIdx,rIdx)]->sigT(eIndx);

	// eventually there will need to be some manipulation here that 
	// extrapolates the cross section based on temperature

	return sigT;
};

//==============================================================================

//==============================================================================
//! sigS return scattering cross section

// return scatter cross section at location indicated by input indices
double Materials::sigS(int zIdx,int rIdx,int gprime,int g){

	double sigS = matBank[matMap(zIdx,rIdx)]->sigS(gprime,g);

	// eventually there will need to be some manipulation here that 
	// extrapolates the cross section based on temperature

	return sigS;
};

//==============================================================================

//==============================================================================
//! sigF return fission cross section 

// return fission cross section at location indicated by input indices
double Materials::sigF(int zIdx,int rIdx,int eIndx){

	double sigF = matBank[matMap(zIdx,rIdx)]->sigF(eIndx);

	// eventually there will need to be some manipulation here that 
	// extrapolates the cross section based on temperature

	return sigF;
};

//==============================================================================

//==============================================================================
//! chiP return probability that a prompt neutron will be born

// return probability at location and energy indicated by input indices
double Materials::chiP(int zIdx,int rIdx,int eIndx){

	double chiP = matBank[matMap(zIdx,rIdx)]->chiP(eIndx);

	// eventually there will need to be some manipulation here that 
	// extrapolates the cross section based on temperature

	return chiP;
};

//==============================================================================

//==============================================================================
//! chiD return probability that a delayed neutron will be born

// return probability at location and energy indicated by input indices
double Materials::chiD(int zIdx,int rIdx,int eIndx){

	double chiD = matBank[matMap(zIdx,rIdx)]->chiD(eIndx);

	// eventually there will need to be some manipulation here that 
	// extrapolates the cross section based on temperature

	return chiD;
};

//==============================================================================

//==============================================================================
//! nu return nu, the average number of neutrons produced per fission event

// return nu at location indicated by input indices
double Materials::nu(int zIdx,int rIdx){

	double nu = matBank[matMap(zIdx,rIdx)]->nu;

	// eventually there will need to be some manipulation here that 
	// extrapolates the cross section based on temperature

	return nu;
};

//==============================================================================

//==============================================================================
//! checkMats Check to see if nuclear data defined on each group is sensible

void Materials::checkMats()
{
	for (int iCount = 0; iCount < matBank.size(); ++iCount){
		matBank[iCount]->checkMat();
	}
};

//==============================================================================

//==============================================================================
//! edit Print out material map and all entries in material bank

void Materials::edit()
{
	cout << "Material map: " << endl;
	cout << matMap << endl;
	cout << endl;
	for (int iCount = 0; iCount < matBank.size(); ++iCount){
		matBank[iCount]->edit();
	}
};

//==============================================================================



