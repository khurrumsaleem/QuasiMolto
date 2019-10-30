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

using namespace std; 
using namespace arma;

//==============================================================================
//! Material class that holds material and geometry information

class Materials
{
        public:
	umat matMap;
        // public functions
        Materials(Mesh * myMesh,YAML::Node * myInput);
	void readMats();
	void readGeom();
        void setMatRegion(int myIndex,double rIn,double rOut,\
		double zUp,double zLow);

        private:
        // private functions
        YAML::Node * input;
        Mesh * mesh;
	map<string,int> mat2idx;
};

//==============================================================================

//==============================================================================
//! Material class object constructor

Materials::Materials(Mesh * myMesh,\
                             YAML::Node * myInput)
{
	// Point to variables for mesh and input file
	mesh = myMesh;
	input = myInput;
	matMap.set_size(mesh->zCent.size(),mesh->rCent.size());
	matMap.fill(0);
	readMats();
	readGeom();
};

//==============================================================================

//==============================================================================
//! readInput Parse input and define  geometry parameters

void Materials::readGeom()
{
	YAML::Node geometry = (*input)["geometry"];
	string regType;
	string matName;
	double rIn, rOut, zLow, zUp;
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
	string regType;
        int iCount=0;
	for (YAML::const_iterator it=mats.begin(); it!=mats.end(); ++it){
		mat2idx.insert(pair<string,int>(it->first.as<string>(),iCount));
		++iCount;
	}
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

