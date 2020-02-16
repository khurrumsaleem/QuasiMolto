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
/// Calculate Eddington factors using angular fluxes from transport objects
void TransportToQDCoupling::calcBCs()
{
  int rows = MGT->SGTs[0]->sFlux.rows();
  int cols = MGT->SGTs[0]->sFlux.cols();
  int eIdx = rows - 1,sIdx = cols - 1;
  int angIdx,xiIdx=0,muIdx=1,etaIdx=2,weightIdx = 3;  
  double angFlux,angFluxN,angFluxS,mu,xi,weight;
  double Jr,JzN,JzS;
  
  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
      for (int iZ = 0; iZ < cols; iZ++)
      {

        // reset accumulators
        Jr = 0.0;

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
            
            Jr = Jr + mu*angFlux*weight;
          
          } //iMu
        } //iXi 

        // set radial BCs 
        MGQD->SGQDs[iGroup]->eCurrentRBC(iZ) = Jr;
        MGQD->SGQDs[iGroup]->wCurrentRBC(iZ) = 0.0;
        
        MGQD->SGQDs[iGroup]->eFluxBC(iZ) \
          = MGT->SGTs[iGroup]->sFlux(iZ,eIdx);
        MGQD->SGQDs[iGroup]->wFluxBC(iZ) \
          = MGT->SGTs[iGroup]->sFlux(iZ,0);

      } //iZ
      for (int iR = 0; iR < cols; iR++)
      {
        // reset accumulators
        JzN = 0.0;
        JzS = 0.0;

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
            
            JzN = JzN + xi*angFluxN*weight;
            JzS = JzS + xi*angFluxS*weight;
          
          } //iMu
        } //iXi 

        // set axial BCs 
        MGQD->SGQDs[iGroup]->nCurrentZBC(iR) = JzN;
        MGQD->SGQDs[iGroup]->sCurrentZBC(iR) = JzS;
        
        MGQD->SGQDs[iGroup]->nFluxBC(iR) \
          = MGT->SGTs[iGroup]->sFlux(0,iR);
        MGQD->SGQDs[iGroup]->sFluxBC(iR) \
          = MGT->SGTs[iGroup]->sFlux(sIdx,iR);

      } //iR
  } //iGroup
}
//==============================================================================