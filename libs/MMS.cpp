#include "MMS.h"

using namespace std;

MMS::MMS(
  MultiGroupTransport * myMGT,\
  Mesh * myMesh,\
  Materials * myMaterials,\
  YAML::Node * myInput)
{
  MGT = myMGT;
  mesh = myMesh;
  input = myInput;
  materials = myMaterials;
  
  if ((*input)["parameters"]["c"]){ 
    c=(*input)["parameters"]["c"].as<double>();
  }
}

Eigen::MatrixXd MMS::isotropicTransportSourceMMS(double xi, double mu, double t) {

  Eigen::MatrixXd mmsSource = Eigen::MatrixXd::Zero(mesh->dzsCorner.size(),\
    mesh->drsCorner.size());

  double zLim = mesh->Z,rLim = mesh->R;
  double zDown,zUp,rDown,rUp,vol;
  double rad1,rad2,rad3,rad4;
  double term1,term2,term3,term4;
  double sinUp,sinDown,cosUp,cosDown;

  double A;
  double B = exp(c*t) * mu;
  double D = exp(c*t) * M_PI * xi/ zLim;

  for (int iZ = 0; iZ < mmsSource.rows(); ++iZ){
    for (int iR = 0; iR < mmsSource.cols(); ++iR){
  
      A = exp(c*t) * (materials->sigT(iZ,iR,0)\
      + c/materials->neutV(0) -1.0*(materials->sigS(iZ,iR,0,0)\
      + materials->nu(iZ,iR,0) * materials->sigF(iZ,iR,0)));
      
      zDown = mesh->zCornerEdge(iZ);
      zUp = mesh->zCornerEdge(iZ+1);
      rDown = mesh->rCornerEdge(iR);
      rUp = mesh->rCornerEdge(iR+1);

      sinUp = sin(M_PI*zUp/zLim); 
      sinDown = sin(M_PI*zDown/zLim); 
      cosUp = cos(M_PI*zUp/zLim); 
      cosDown = cos(M_PI*zDown/zLim); 

      rad1 = ((pow(rLim,2)*pow(rUp,2)/2.0-pow(rUp,4)/4.0)
        -(pow(rLim,2)*pow(rDown,2)/2.0 - pow(rDown,4)/4.0));
      term1 = -(A*zLim/M_PI)*(cosUp-cosDown)*rad1;

      rad2 = ((pow(rLim,2)*rUp-pow(rUp,3)/3.0)\
        -(pow(rLim,2)*rDown-pow(rDown,3)/3.0));
      term2 = (B*zLim/M_PI)*(cosUp-cosDown)*rad2;

      rad3 = ((pow(rLim,2)*rUp-pow(rUp,3))\
        -(pow(rLim,2)*rDown - pow(rDown,3)));
      term3 = -(B*zLim/M_PI)*(cosUp-cosDown)*rad3;

      rad4 = ((pow(rLim,2)*pow(rUp,2)/2.0-pow(rUp,4)/4.0)\
        -(pow(rLim,2)*pow(rDown,2)/2.0-pow(rDown,4)/4.0));
      term4 = (D*zLim/M_PI)*(sinUp-sinDown)*rad4;
     
      vol = (zUp - zDown) * (pow(rUp,2)-pow(rDown,2))/2.0; 

      mmsSource(iZ,iR) = (term1 + term2 + term3 + term4)/vol;      
    }
  }

  return mmsSource;
}

Eigen::MatrixXd MMS::anisotropicTransportSourceMMS(double xi, double mu, double t) {

  Eigen::MatrixXd mmsSource = Eigen::MatrixXd::Zero(mesh->dzsCorner.size(),\
    mesh->drsCorner.size());

  double eta,etaSquared = 1.0-pow(xi,2)-pow(mu,2);
  
  if (etaSquared < 1E-2){
    eta = 0.0;
  } else {
    eta = sqrt(etaSquared);
  } 

  double zLim = mesh->Z,rLim = mesh->R;
  double zDown,zUp,rDown,rUp,vol;
  double rad1,rad2,rad3,rad4;
  double term1,term2,term3,term4,term5;
  double sinUp,sinDown,cosUp,cosDown;

  double A = exp(c*t) * pow(mu,3);
  double B = exp(c*t) * M_PI * xi * pow(mu,2) / zLim;
  double D = -exp(c*t) * (pow(mu,3)-2*pow(eta,2)*mu);
  double F,relUp,relDown;

  for (int iZ = 0; iZ < mmsSource.rows(); ++iZ){
    for (int iR = 0; iR < mmsSource.cols(); ++iR){
      
      F = exp(c*t) * ( pow(mu,2) * (materials->sigT(iZ,iR,0)\
      + c/materials->neutV(0)) - (1.0/3.0)\
      *(materials->sigS(iZ,iR,0,0)
      +materials->nu(iZ,iR,0) * materials->sigF(iZ,iR,0)));
  
      zDown = mesh->zCornerEdge(iZ);
      zUp = mesh->zCornerEdge(iZ+1);
      rDown = mesh->rCornerEdge(iR);
      rUp = mesh->rCornerEdge(iR+1);

      sinUp = sin(M_PI*zUp/zLim); 
      sinDown = sin(M_PI*zDown/zLim); 
      cosUp = cos(M_PI*zUp/zLim); 
      cosDown = cos(M_PI*zDown/zLim); 

      relUp = 2.0*M_PI*rUp/rLim;
      relDown = 2.0*M_PI*rDown/rLim;

      rad1 = rUp*(1.0-cos(relUp))-rDown*(1.0-cos(relDown));
      term1 = -(A*zLim/M_PI)*(cosUp-cosDown)*rad1;

      rad2 = pow(rUp,2.0)/2.0\
      - rLim*(rLim*cos(relUp)+2.0*M_PI*rUp*sin(relUp))/(4.0*pow(M_PI,2.0))\
      - (pow(rDown,2.0)/2.0\
      - rLim*(rLim*cos(relDown)+2.0*M_PI*rDown*sin(relDown))/(4.0*pow(M_PI,2.0)));
      term2 = (B*zLim/M_PI)*(sinUp-sinDown)*rad2;

      rad3 = rUp - rLim*sin(relUp)/(2.0*M_PI)\
      - (rDown - rLim*sin(relDown)/(2.0*M_PI));
      term3 = -(D*zLim/M_PI)*(cosUp-cosDown)*rad3;

      term4 = -(F*zLim/M_PI)*(cosUp-cosDown)*rad2;

      vol = (zUp - zDown) * (pow(rUp,2.0)-pow(rDown,2))/2.0; 
      
      mmsSource(iZ,iR) = (term1 + term2 + term3 + term4)/vol;      
    }
  }

  return mmsSource;
}


void MMS::timeDependent(){

  Eigen::MatrixXd isoSource;
  double xi,mu;
  double t = mesh->dt;
  double fluxResidual=0,alphaResidual=0,fissionResidual=0;

  // set first set of alphas
  MGT->SGTs[0]->alpha.setOnes();  
  MGT->SGTs[0]->alpha = c*MGT->SGTs[0]->alpha;  
  
  // compute source isotropic
  MGT->SGTs[0]->calcSource();

  isoSource = MGT->SGTs[0]->q;
  
  // loop over angles 
  for (int idt = 0; idt < 1; ++idt){
    t = idt*mesh->dt+mesh->dt;    

    for (int sourceIter = 0; sourceIter < 100; ++sourceIter){

      for (int scatterIter = 0; scatterIter < 1000; ++scatterIter){
        
        // first solve for all starting angles 
        for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi){
          xi = mesh->quadrature[iXi].quad[0][0];
          mu = -sin(acos(xi));
          MGT->SGTs[0]->q = isoSource + anisotropicTransportSourceMMS(xi,mu,t);

          // call solver for each angle
          MGT->startAngleSolve->solveAngularFlux(&(MGT->SGTs[0]->aHalfFlux),\
          &(MGT->SGTs[0]->q),\
          &(MGT->SGTs[0]->alpha),\
          0,\
          iXi);

        }

        for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi){
          xi = mesh->quadrature[iXi].quad[0][0];
            
          for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu){
            mu = mesh->quadrature[iXi].quad[iMu][1];
     
            // calculate source at t1 for each angle
            MGT->SGTs[0]->q = isoSource + anisotropicTransportSourceMMS(xi,mu,t);
            
            // call solver for each angle
            MGT->SCBSolve->solveAngularFlux(&(MGT->SGTs[0]->aFlux),\
            &(MGT->SGTs[0]->aHalfFlux),\
            &(MGT->SGTs[0]->q),\
            &(MGT->SGTs[0]->alpha),\
            0,\
            iXi,\
            iMu);
          } //iMu
        } //iXi
    
        // calculate scalar flux
        fluxResidual = MGT->SGTs[0]->calcFlux();
        cout << "Flux residual: " << endl;
        cout << fluxResidual << endl;

        // calculate scatter source
        MGT->SGTs[0]->calcSource("s");
  
        // store new scattering source
        isoSource = MGT->SGTs[0]->q;
        
        if (fluxResidual < MGT->epsFlux)
          break;
      }
      
      fissionResidual = MGT->SGTs[0]->calcSource();
        
      // store new source
      isoSource = MGT->SGTs[0]->q;
        
      cout << "Fission residual: " << endl;
      cout << fissionResidual << endl;

      if (fissionResidual < MGT->epsFissionSource){
        cout << "Fission convergence:" << MGT->epsFissionSource << endl;
        if (idt != 0){
          alphaResidual = MGT->SGTs[0]->calcAlpha();
          cout << "Alpha residual: " << endl;
          cout << alphaResidual << endl;
          if (alphaResidual < MGT->epsAlpha)
            break;
        }
        else
          break;
      }
    }
 
    // set sFluxPrev
    MGT->SGTs[0]->sFluxPrev = MGT->SGTs[0]->sFlux;
  
    // set first set of alphas
    MGT->SGTs[0]->alpha.setOnes();  
    MGT->SGTs[0]->alpha = c*MGT->SGTs[0]->alpha;  
     
  } //idt

  MGT->SGTs[0]->writeFlux();

}

