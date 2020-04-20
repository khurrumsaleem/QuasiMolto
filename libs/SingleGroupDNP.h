#ifndef SingleGroupDNP_H
#define SingleGroupDNP_H

#include "Mesh.h"
#include "Materials.h"

using namespace std;

class MultiGroupDNP; // forward declaration

//==============================================================================
//! Contains information and builds linear system for a single precursor group

class SingleGroupDNP
{
  public:
  Eigen::MatrixXd dnpConc;
  double beta; 
  double lambda;
  int indexOffset = 0;
  SingleGroupDNP(Materials * myMats,\
    Mesh * myMesh,\
    MultiGroupDNP * myMGDNPS,\
    double myBeta,\
    double myLambda);
  int getIndex(int iZ,int iR);
  void calcDiracs();
  void calcFluxes();
  void getConc();
  double calcPhi(double theta,string fluxLimiter); 
  double calcTheta(double TupwindInterface,double Tinterface);

  private: 
  Materials * mats;
  Mesh * mesh;
  MultiGroupDNP * mgdnp;
};

//==============================================================================

#endif
