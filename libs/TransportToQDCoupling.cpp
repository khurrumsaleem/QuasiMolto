// File: TransportToQDCoupling.cpp
// Purpose: handle communication between transport and quasidiffusion objects
// Date: February 12, 2020

#include <iostream>
#include <vector>
#include <iomanip>
#include <armadillo>
#include "Mesh.h"
#include "Materials.h"
#include "MultiGroupQD.h"
#include "MultiGroupTransport.h"
#include "SingleGroupQD.h"
#include "SingleGroupTransport.h"
#include "TransportToQDCoupling.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../TPLs/eigen-git-mirror/Eigen/Eigen"

using namespace std;
using namespace arma;

//==============================================================================
/// TransportToQDCoupling class object constructor
///
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myInput YAML input object for the simulation
TransportToQDCoupling::TransportToQDCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput,\
  MultiGroupTransport * myMGT,\
  MultiGroupQD * myMGQD)
{
  // Assign pointers for materials, mesh, and input objects
  materials = myMaterials;
  mesh = myMesh;
  input = myInput;

  // Assign pointers for multigroup transport and quasidiffusion objects
  MGT = myMGT;
  MGQD = myMGQD;

};
//==============================================================================

//==============================================================================
/// Calculate Eddington factors using angular fluxes from transport objects
void TransportToQDCoupling::calcEddingtonFactors()
{
  int rows = MGT->SGTs[0]->sFlux.rows();
  int cols = MGT->SGTs[0]->sFlux.cols();
  int angIdx,xiIdx=0,muIdx=1,etaIdx=2,weightIdx = 3;  
  double angFlux,mu,xi,weight,EzzCoef,ErrCoef,ErzCoef;
  double numeratorEzz,numeratorErr,numeratorErz,denominator; 
  
  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    for (int iR = 0; iR < rows; iR++)
    {
      for (int iZ = 0; iZ < cols; iZ++)
      {

        // reset accumulators
        numeratorEzz = 0.0;
        numeratorErr = 0.0;
        numeratorErz = 0.0;
        denominator = 0.0;

        // loop over quadrature
        for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi)
        {
          xi = mesh->quadrature[iXi].quad[0][xiIdx];
          for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu)
          {
            angIdx=mesh->quadrature[iXi].ordIdx[iMu];
            mu = mesh->quadrature[iXi].quad[iMu][muIdx]; 
            weight = mesh->quadrature[iXi].quad[iMu][weightIdx]; 
            angFlux = MGT->SGTs[iGroup]->aFlux(iZ,iR,angIdx);

            EzzCoef = xi*xi;
            ErrCoef = mu*mu;
            ErzCoef = mu*xi;

            numeratorEzz = numeratorEzz + EzzCoef*angFlux*weight;
            numeratorErr = numeratorErr + ErrCoef*angFlux*weight;
            numeratorErz = numeratorErz + ErzCoef*angFlux*weight;

            denominator = denominator + angFlux*weight;
          
          } //iMu
        } //iXi 
 
        MGQD->SGQDs[iGroup]->Ezz(iZ,iR) = numeratorEzz/denominator;
        MGQD->SGQDs[iGroup]->Err(iZ,iR) = numeratorErr/denominator;
        MGQD->SGQDs[iGroup]->Erz(iZ,iR) = numeratorErz/denominator;

      } //iZ
    } //iR
    cout << "iGroup: " << iGroup << endl;
    cout << "Ezz: " << endl;
    cout << MGQD->SGQDs[iGroup]->Ezz << endl;
    cout << "Err: " << endl;
    cout << MGQD->SGQDs[iGroup]->Err << endl;
    cout << "Erz: " << endl;
    cout << MGQD->SGQDs[iGroup]->Erz << endl;
  } //iGroup
}
//==============================================================================

//==============================================================================
/// Calculate boundary conditions for scalar flux 
void TransportToQDCoupling::calcBCs()
{
  int rows = MGT->SGTs[0]->sFlux.rows();
  int cols = MGT->SGTs[0]->sFlux.cols();
  int eIdx = rows - 1,sIdx = cols - 1;
  int angIdx,xiIdx=0,muIdx=1,etaIdx=2,weightIdx = 3;  
  double angFlux,angFluxN,angFluxS,mu,xi,weight;
  double inwardJrE,inwardJzN,inwardJzS;
  double inwardFluxE,inwardFluxN,inwardFluxS;
  double outwardJrE,outwardJzN,outwardJzS;
  double outwardFluxE,outwardFluxN,outwardFluxS;
  
  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
      for (int iZ = 0; iZ < cols; iZ++)
      {

        // reset accumulators
        inwardJrE = 0.0;
        inwardFluxE = 0.0;
        outwardJrE = 0.0;
        outwardFluxE = 0.0;

        // loop over quadrature
        for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi)
        {
          xi = mesh->quadrature[iXi].quad[0][xiIdx];
          for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu)
          {
            angIdx=mesh->quadrature[iXi].ordIdx[iMu];
            mu = mesh->quadrature[iXi].quad[iMu][muIdx]; 
            weight = mesh->quadrature[iXi].quad[iMu][weightIdx];
 
            angFlux = MGT->SGTs[iGroup]->aFlux(iZ,eIdx,angIdx);

            // only accumulate inward facing angular fluxes on the the 
            // outside radial (east) boundary 
            if (mu < 0) 
            {
              inwardJrE = inwardJrE + mu*angFlux*weight;
              inwardFluxE = inwardFluxE + angFlux*weight;
            } else
            {
              outwardJrE = outwardJrE + mu*angFlux*weight;
              outwardFluxE = outwardFluxE + angFlux*weight;
            }
          } //iMu
        } //iXi 

        // set inward current in SGQD object 
        MGQD->SGQDs[iGroup]->eInwardCurrentBC(iZ) = inwardJrE;

        // set inward flux in SGQD object 
        MGQD->SGQDs[iGroup]->eInwardFluxBC(iZ) = inwardFluxE;

        // set outward current to flux ratio in SGQD object 
        MGQD->SGQDs[iGroup]->eOutwardCurrToFluxRatioBC(iZ)\
          = outwardJrE/outwardFluxE;

      } //iZ
      for (int iR = 0; iR < cols; iR++)
      {
        // reset accumulators
        inwardJzN = 0.0;
        inwardJzS = 0.0;
        inwardFluxN = 0.0;
        inwardFluxS = 0.0;
        outwardJzN = 0.0;
        outwardJzS = 0.0;
        outwardFluxN = 0.0;
        outwardFluxS = 0.0;

        // loop over quadrature
        for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi)
        {
          xi = mesh->quadrature[iXi].quad[0][xiIdx];
          for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu)
          {
            angIdx=mesh->quadrature[iXi].ordIdx[iMu];
            mu = mesh->quadrature[iXi].quad[iMu][muIdx]; 
            weight = mesh->quadrature[iXi].quad[iMu][weightIdx];
 
            angFluxN = MGT->SGTs[iGroup]->aFlux(0,iR,angIdx);
            angFluxS = MGT->SGTs[iGroup]->aFlux(sIdx,iR,angIdx);
  
            // accumulate outward and inward angular fluxes in separate
            // variables on the north face
            if (xi > 0) 
            {
              inwardJzN = inwardJzN + xi*angFluxN*weight;
              inwardFluxN = inwardFluxN + angFluxN*weight;
            } else
            {
              outwardJzN = outwardJzN + xi*angFluxN*weight;
              outwardFluxN = outwardFluxN + angFluxN*weight;
            }
            // accumulate outward and inward angular fluxes in separate
            // variables on the south face
            if (xi < 0)
            {
              inwardJzS = inwardJzS + xi*angFluxS*weight;
              inwardFluxS = inwardFluxS + angFluxS*weight;
            } else
            {
              outwardJzS = outwardJzS + xi*angFluxS*weight;
              outwardFluxS = outwardFluxS + angFluxS*weight;
            }
          
          } //iMu
        } //iXi 

        // set inward current in SGQD object 
        MGQD->SGQDs[iGroup]->nInwardCurrentBC(iR) = inwardJzN;
        MGQD->SGQDs[iGroup]->sInwardCurrentBC(iR) = inwardJzS;
        
        // set inward flux in SGQD object 
        MGQD->SGQDs[iGroup]->nInwardFluxBC(iR) = inwardFluxN;
        MGQD->SGQDs[iGroup]->sInwardFluxBC(iR) = inwardFluxS;
  
        // set outward current to flux ratio in SGQD object
        MGQD->SGQDs[iGroup]->nOutwardCurrToFluxRatioBC(iR)\
          = outwardJzN/outwardFluxN;
        MGQD->SGQDs[iGroup]->sOutwardCurrToFluxRatioBC(iR)\
          = outwardJzS/outwardFluxS;

      } //iR  
  
  } //iGroup

}
//==============================================================================
