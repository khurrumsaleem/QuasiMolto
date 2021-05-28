// File: HeatTransfer.cpp     
// Purpose: Contains information and builds linear system for heat transfer
// Date: April 9, 2020

#include "HeatTransfer.h"
#include "MultiPhysicsCoupledQD.h"

using namespace std;

//==============================================================================
/// HeatTransfer class object constructor
///
/// @param [in] myMaterials materials object for the simulation
/// @param [in] myMesh mesh object for the simulation
/// @param [in] myInput input object for the simulation
/// @param [in] myQD multiphysics coupled QD object for the simulation
HeatTransfer::HeatTransfer(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * myInput,\
  MultiPhysicsCoupledQD * myMPQD)
{
  // Assign inputs to their member variables
  mats = myMaterials;
  mesh = myMesh;
  input = myInput;
  mpqd = myMPQD;

  checkOptionalParams();

  // Set number of unknowns
  nUnknowns = mesh->nR*mesh->nZ;

  // Initialize size of matrices and vectors
  temp.setConstant(mesh->nZ,mesh->nR,inletT);
  flux.setZero(mesh->nZ+1,mesh->nR);
  dirac.setZero(mesh->nZ+1,mesh->nR);
  inletTemp.setConstant(2,mesh->nR,inletT);
  outletTemp.setZero(mesh->nR);
  inletDensity.setZero(mesh->nR);
  inletVelocity.setZero(mesh->nR);
  inletcP.setZero(mesh->nR);
  
  assignBoundaryIndices();
};
//==============================================================================

//==============================================================================
/// Build linear system to solve heat equation
///
void HeatTransfer::buildLinearSystem()
{
  
  int myIndex,sIndex,nIndex,wIndex,eIndex,iEq = indexOffset;
  int iEqTemp = 0;
  int nR = temp.cols()-1;
  int nZ = temp.rows()-1;
  double harmonicAvg,coeff,cCoeff;
  Eigen::MatrixXd volAvgGammaDep;
  vector<double> gParams;

  updateBoundaryConditions();
  calcDiracs();
  calcFluxes();

  // Calculate core-average gamma deposition term
  if (modIrradiation == axial)
    volAvgGammaDep = calcExplicitAxialFissionEnergy();
  else if (modIrradiation == fuel)
    volAvgGammaDep = calcExplicitAxialFuelFissionEnergy();
  else
    volAvgGammaDep = calcExplicitFissionEnergy();

  
  //Atemp.resize(nUnknowns,mpqd->A.cols());
  //Atemp.reserve(5*nUnknowns);
  
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Atemp;
  Atemp.resize(nUnknowns,mpqd->A.cols());
  Atemp.setZero();
  
  #pragma omp parallel for private(myIndex,sIndex,nIndex,wIndex,eIndex,\
    gParams,cCoeff,coeff,harmonicAvg,iEq,iEqTemp)
  for (int iZ = 0; iZ < temp.rows(); iZ++)
  {

    for (int iR = 0; iR < temp.cols(); iR++)
    {
      
      iEq = getIndex(iZ,iR);
      iEqTemp = iEq - indexOffset;

      // Reset center coefficient
      cCoeff = 0;

      // Get cell indices
      myIndex = getIndex(iZ,iR);     
      sIndex = getIndex(iZ+1,iR);     
      nIndex = getIndex(iZ-1,iR);     
      wIndex = getIndex(iZ,iR-1);     
      eIndex = getIndex(iZ,iR+1);     

      gParams = mesh->getGeoParams(iR,iZ);
      
      //Atemp.insert(iEqTemp,myIndex) = mats->density(iZ,iR)*mats->cP(iZ,iR);
      cCoeff = mats->density(iZ,iR)*mats->cP(iZ,iR);
      
      // East face
      if (iR == nR)
      {
        coeff = (-gParams[iEF]*mats->k(iZ,iR)\
        /mesh->drsCorner(iR))/gParams[iVol];
        //Atemp.coeffRef(iEqTemp,myIndex) -= mesh->dt*coeff;
        cCoeff -= mesh->dt*coeff;
        mpqd->b(iEq) -= mesh->dt*coeff*wallT; 
      } else
      {
        harmonicAvg = pow(mesh->drsCorner(iR)/mats->k(iZ,iR)\
          + mesh->drsCorner(iR+1)/mats->k(iZ,iR+1),-1.0);
        coeff = -2.0*gParams[iEF]*harmonicAvg/gParams[iVol];
        //Atemp.insert(iEqTemp,eIndex) = mesh->dt*coeff;
        Atemp(iEqTemp,eIndex) = mesh->dt*coeff;
        //Atemp.coeffRef(iEqTemp,myIndex) -= mesh->dt*coeff;
        cCoeff -= mesh->dt*coeff;
      }

      // West face
      if (iR != 0)
      {
        harmonicAvg = pow(mesh->drsCorner(iR-1)/mats->k(iZ,iR-1)\
          + mesh->drsCorner(iR)/mats->k(iZ,iR),-1.0);
        coeff = 2.0*gParams[iWF]*harmonicAvg/gParams[iVol];
        //Atemp.insert(iEqTemp,wIndex) = -mesh->dt*coeff;
        Atemp(iEqTemp,wIndex) = -mesh->dt*coeff;
        //Atemp.coeffRef(iEqTemp,myIndex) += mesh->dt*coeff;
        cCoeff += mesh->dt*coeff;
      } 

      // North face
      if (iZ == 0 and mats->posVelocity)
      {
        coeff = (gParams[iNF]*mats->k(iZ,iR)\
        /mesh->dzsCorner(iZ))/gParams[iVol];
        //Atemp.coeffRef(iEqTemp,myIndex) += mesh->dt*coeff;
        cCoeff += mesh->dt*coeff;
        mpqd->b(iEq) += mesh->dt*coeff*inletTemp(1,iR);             
      } else if (iZ != 0)
      {
        harmonicAvg = pow(mesh->dzsCorner(iZ-1)/mats->k(iZ-1,iR)\
          + mesh->dzsCorner(iZ)/mats->k(iZ,iR),-1.0);
        coeff = 2.0*gParams[iNF]*harmonicAvg/gParams[iVol];
        //Atemp.insert(iEqTemp,nIndex) = -mesh->dt*coeff;
        Atemp(iEqTemp,nIndex) = -mesh->dt*coeff;
        //Atemp.coeffRef(iEqTemp,myIndex) += mesh->dt*coeff;
        cCoeff += mesh->dt*coeff;
      }

      // South face
      if (iZ == nZ and !(mats->posVelocity))
      {
        coeff = -(gParams[iSF]*mats->k(iZ,iR)\
        /mesh->dzsCorner(iZ))/gParams[iVol];
        //Atemp.coeffRef(iEqTemp,myIndex) -= mesh->dt*coeff;
        cCoeff -= mesh->dt*coeff;
        mpqd->b(iEq) -= mesh->dt*coeff*inletTemp(0,iR);             
      } else if (iZ != nZ)
      {
        harmonicAvg = pow(mesh->dzsCorner(iZ+1)/mats->k(iZ+1,iR)\
          + mesh->dzsCorner(iZ)/mats->k(iZ,iR),-1.0);
        coeff = -2.0*gParams[iSF]*harmonicAvg/gParams[iVol];
        //Atemp.insert(iEqTemp,sIndex) = mesh->dt*coeff;
        Atemp(iEqTemp,sIndex) = mesh->dt*coeff;
        //Atemp.coeffRef(iEqTemp,myIndex) -= mesh->dt*coeff;
        cCoeff -= mesh->dt*coeff;
      }

      // Insert cell center coefficient
      //Atemp.insert(iEqTemp,myIndex) = cCoeff;
      Atemp(iEqTemp,myIndex) = cCoeff;

      // Time term
      mpqd->b(iEq) += mats->density(iZ,iR)*mats->cP(iZ,iR)*temp(iZ,iR); 

      // Flux source term 
      coeff = -mesh->dt*mats->omega(iZ,iR)*mats->oneGroupXS->sigF(iZ,iR);
      //mpqd->fluxSource(iZ,iR,iEqTemp,coeff,&Atemp);
      mpqd->fluxSource(iZ,iR,iEqTemp,coeff,&Atemp);
      
      // Gamma source term 
      coeff = mesh->dt;
      mpqd->b(iEq) += coeff * mats->gamma(iZ,iR) * volAvgGammaDep(iZ,iR);
      //gammaSource(iZ,iR,iEqTemp,coeff);
      

      // Advection term
      mpqd->b(iEq) += (mesh->dt/mesh->dzsCorner(iZ))*(flux(iZ,iR)-flux(iZ+1,iR));

      // Iterate equation count
      //iEq = iEq + 1;
      //iEqTemp = iEqTemp + 1;
  
    }
  }
   
  mpqd->A.middleRows(indexOffset,nUnknowns) = Atemp.sparseView(); 

};
//==============================================================================

//==============================================================================
/// Build steady state linear system to solve heat equation
///
void HeatTransfer::buildSteadyStateLinearSystem()
{
  
  int myIndex,sIndex,nIndex,wIndex,eIndex,upwindIndex,iEq = indexOffset;
  int iEqTemp = 0;
  int nR = temp.cols()-1;
  int nZ = temp.rows()-1;
  double harmonicAvg,coeff,keff,neutronFlux,cCoeff;
  Eigen::MatrixXd volAvgGammaDep;
  vector<double> gParams;

  updateBoundaryConditions();
  calcImplicitFluxes();

  // Calculate core-average gamma deposition term
  if (modIrradiation == axial)
    volAvgGammaDep = calcExplicitAxialFissionEnergy();
  else if (modIrradiation == fuel)
    volAvgGammaDep = calcExplicitAxialFuelFissionEnergy();
  else
    volAvgGammaDep = calcExplicitFissionEnergy();
  
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Atemp;
  Atemp.resize(nUnknowns,mpqd->A.cols());
  Atemp.setZero();
  
  #pragma omp parallel for private(myIndex,sIndex,nIndex,wIndex,eIndex,\
    upwindIndex,gParams,cCoeff,coeff,keff,neutronFlux,harmonicAvg,iEq,iEqTemp)
  for (int iZ = 0; iZ < temp.rows(); iZ++)
  {

    for (int iR = 0; iR < temp.cols(); iR++)
    {
      
      iEq = getIndex(iZ,iR);
      iEqTemp = iEq - indexOffset;

      // Reset center coefficient
      cCoeff = 0;

      // Get cell indices
      myIndex = getIndex(iZ,iR);     
      sIndex = getIndex(iZ+1,iR);     
      nIndex = getIndex(iZ-1,iR);     
      wIndex = getIndex(iZ,iR-1);     
      eIndex = getIndex(iZ,iR+1);     

      gParams = mesh->getGeoParams(iR,iZ);
      
      // East face
      if (iR == nR)
      {
        coeff = (-gParams[iEF]*mats->k(iZ,iR)\
        /mesh->drsCorner(iR))/gParams[iVol];
        //cCoeff -= mesh->dt*coeff;
        cCoeff -= coeff;
        //mpqd->b(iEq) -= mesh->dt*coeff*wallT; 
        mpqd->b(iEq) -= coeff*wallT; 
      } else
      {
        harmonicAvg = pow(mesh->drsCorner(iR)/mats->k(iZ,iR)\
          + mesh->drsCorner(iR+1)/mats->k(iZ,iR+1),-1.0);
        coeff = -2.0*gParams[iEF]*harmonicAvg/gParams[iVol];
        //Atemp(iEqTemp,eIndex) = mesh->dt*coeff;
        Atemp(iEqTemp,eIndex) = coeff;
        //cCoeff -= mesh->dt*coeff;
        cCoeff -= coeff;
      }

      // West face
      if (iR != 0)
      {
        harmonicAvg = pow(mesh->drsCorner(iR-1)/mats->k(iZ,iR-1)\
          + mesh->drsCorner(iR)/mats->k(iZ,iR),-1.0);
        coeff = 2.0*gParams[iWF]*harmonicAvg/gParams[iVol];
        //Atemp(iEqTemp,wIndex) = -mesh->dt*coeff;
        Atemp(iEqTemp,wIndex) = -coeff;
        //cCoeff += mesh->dt*coeff;
        cCoeff += coeff;
      } 

      // North face
      if (iZ == 0 and mats->posVelocity)
      {
        coeff = (gParams[iNF]*mats->k(iZ,iR)\
        /mesh->dzsCorner(iZ))/gParams[iVol];
        //cCoeff += mesh->dt*coeff;
        cCoeff += coeff;
        //mpqd->b(iEq) += mesh->dt*coeff*inletTemp(1,iR);             
        mpqd->b(iEq) += coeff*inletTemp(1,iR);             
      } else if (iZ != 0)
      {
        harmonicAvg = pow(mesh->dzsCorner(iZ-1)/mats->k(iZ-1,iR)\
          + mesh->dzsCorner(iZ)/mats->k(iZ,iR),-1.0);
        coeff = 2.0*gParams[iNF]*harmonicAvg/gParams[iVol];
        //Atemp(iEqTemp,nIndex) = -mesh->dt*coeff;
        Atemp(iEqTemp,nIndex) = -coeff;
        //cCoeff += mesh->dt*coeff;
        cCoeff += coeff;
      }

      // South face
      if (iZ == nZ and !(mats->posVelocity))
      {
        coeff = -(gParams[iSF]*mats->k(iZ,iR)\
        /mesh->dzsCorner(iZ))/gParams[iVol];
        //cCoeff -= mesh->dt*coeff;
        cCoeff -= coeff;
        //mpqd->b(iEq) -= mesh->dt*coeff*inletTemp(0,iR);             
        mpqd->b(iEq) -= coeff*inletTemp(0,iR);             
      } else if (iZ != nZ)
      {
        harmonicAvg = pow(mesh->dzsCorner(iZ+1)/mats->k(iZ+1,iR)\
          + mesh->dzsCorner(iZ)/mats->k(iZ,iR),-1.0);
        coeff = -2.0*gParams[iSF]*harmonicAvg/gParams[iVol];
        //Atemp(iEqTemp,sIndex) = mesh->dt*coeff;
        Atemp(iEqTemp,sIndex) = coeff;
        //cCoeff -= mesh->dt*coeff;
        cCoeff -= coeff;
      }

      // Insert cell center coefficient
      //Atemp.insert(iEqTemp,myIndex) = cCoeff;
      Atemp(iEqTemp,myIndex) = cCoeff;

      // Time term
      //mpqd->b(iEq) += mats->density(iZ,iR)*mats->cP(iZ,iR)*temp(iZ,iR); 

      // Flux source term 
      coeff = mats->omega(iZ,iR)*mats->oneGroupXS->sigF(iZ,iR);
      keff = mats->oneGroupXS->keff; 
      neutronFlux = mpqd->ggqd->sFlux(iZ,iR);
       
      //mpqd->b(iEq) += coeff*neutronFlux/keff; 
      mpqd->b(iEq) += coeff*neutronFlux; 
      //mpqd->fluxSource(iZ,iR,iEqTemp,coeff,&Atemp);
      //mpqd->fluxSource(iZ,iR,iEqTemp,coeff,&Atemp);
      
      // Gamma source term 
      //coeff = mesh->dt;
      //mpqd->b(iEq) += coeff * mats->gamma(iZ,iR) * volAvgGammaDep;
      mpqd->b(iEq) += mats->gamma(iZ,iR) * volAvgGammaDep(iZ,iR);
      //gammaSource(iZ,iR,iEqTemp,coeff);

      // Advection terms
      //mpqd->b(iEq) += (mesh->dt/mesh->dzsCorner(iZ))*(flux(iZ,iR)-flux(iZ+1,iR));

      if (mats->posVelocity) 
      {
        upwindIndex = getIndex(iZ-1,iR);     

        // Advection terms

        // Upwind cell
        if (iZ == 0) // boundary case
          mpqd->b(iEq) += flux(iZ,iR)*inletTemp(1,iR)/mesh->dzsCorner(iZ);
        else
          Atemp.coeffRef(iEqTemp,upwindIndex) += -flux(iZ,iR)/mesh->dzsCorner(iZ);

        // Primary cell
        Atemp.coeffRef(iEqTemp,myIndex) += flux(iZ+1,iR)/mesh->dzsCorner(iZ);
      }
      else
      {
        upwindIndex = getIndex(iZ+1,iR);     

        // Advection terms

        // Upwind cell
        if (iZ == temp.rows()) // boundary case
          mpqd->b(iEq) -= flux(iZ+1,iR)*inletTemp(0,iR)/mesh->dzsCorner(iZ);
        else
          Atemp.coeffRef(iEqTemp,upwindIndex) = flux(iZ+1,iR)/mesh->dzsCorner(iZ);

        // Primary cell
        Atemp.coeffRef(iEqTemp,myIndex) -= flux(iZ,iR)/mesh->dzsCorner(iZ);
      }
    }
  }

  mpqd->A.middleRows(indexOffset,nUnknowns) = Atemp.sparseView(); 

};
//==============================================================================

//==============================================================================
/// Assert a energy deposition term from gamma rays
///
/// @param [in] iZ axial index 
/// @param [in] iR radial index 
/// @param [in] iEq equation index 
/// @param [in] coeff Coefficient to multiple gamma source term by  
void HeatTransfer::gammaSource(int iZ,int iR,int iEq,double coeff,\
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> * myA)
{

  int myIndex;
  double localVolume,totalVolume;
  double localGamma,localSigF,localOmega,gammaSourceCoeff;
  localGamma = mats->gamma(iZ,iR);

  totalVolume = M_PI*mesh->R*mesh->R*mesh->Z;

  for (int iR = 0; iR < temp.cols(); iR++)
  {
    for (int iZ = 0; iZ < temp.rows(); iZ++)
    {
      // Get local parameters
      localSigF = mats->oneGroupXS->sigF(iZ,iR);
      localOmega = mats->omega(iZ,iR);
      localVolume = mesh->getGeoParams(iR,iZ)[0];
 
      // Calculate gamma source coefficient 
      gammaSourceCoeff = coeff*localGamma*localOmega*localSigF;
      gammaSourceCoeff = localVolume*gammaSourceCoeff/totalVolume; 
      mpqd->fluxSource(iZ,iR,iEq,gammaSourceCoeff,myA);
      
    }
  }
   
};
//==============================================================================

//==============================================================================
/// Calculate core-average gamma energy deposition term
///
Eigen::MatrixXd HeatTransfer::calcExplicitFissionEnergy()
{
  
  double localVolume,localFlux,totalVolume,localSigF,localOmega,\
    volAvgGammaDep = 0,fuelVol = 0;
  
  Eigen::MatrixXd avgFissionEnergy;
 
  totalVolume = M_PI*mesh->R*mesh->R*mesh->Z;
 
  for (int iZ = 0; iZ < temp.rows(); iZ++)
  {
    for (int iR = 0; iR < temp.cols(); iR++)
    {
      // Get local parameters
      localSigF = mats->oneGroupXS->sigF(iZ,iR);
      localOmega = mats->omega(iZ,iR);
      localFlux = mpqd->ggqd->sFlux(iZ,iR);
      localVolume = mesh->getGeoParams(iR,iZ)[0];

      // Calculate gamma source coefficient
      volAvgGammaDep += localVolume*(localOmega*localSigF*localFlux); 
    }
  }

  volAvgGammaDep = volAvgGammaDep/totalVolume; 

  avgFissionEnergy.setConstant(temp.rows(),temp.cols(),volAvgGammaDep);

  return avgFissionEnergy;

};
//==============================================================================

//==============================================================================
/// Calculate core-average gamma energy deposition term
///
Eigen::MatrixXd HeatTransfer::calcExplicitAxialFissionEnergy()
{
  
  double localVolume,localFlux,totalVolume,localSigF,localOmega,\
    axialSliceVolume = 0,axialSliceEnergy = 0,fuelVol = 0;
 
  Eigen::MatrixXd avgFissionEnergy, temporaryMatrix;
  avgFissionEnergy.setZero(temp.rows(),temp.cols());
 
  totalVolume = M_PI*mesh->R*mesh->R*mesh->Z;
 
  for (int iZ = 0; iZ < temp.rows(); iZ++)
  {
    // Reset slice quantities
    axialSliceEnergy = 0;
    axialSliceVolume = 0;

    for (int iR = 0; iR < temp.cols(); iR++)
    {
      // Get local parameters
      localSigF = mats->oneGroupXS->sigF(iZ,iR);
      localOmega = mats->omega(iZ,iR);
      localFlux = mpqd->ggqd->sFlux(iZ,iR);
      localVolume = mesh->getGeoParams(iR,iZ)[0];
      axialSliceVolume += localVolume;

      // Calculate gamma source coefficient
      axialSliceEnergy += localVolume*(localOmega*localSigF*localFlux); 
    }
      
    axialSliceEnergy = axialSliceEnergy/axialSliceVolume;
    temporaryMatrix.setConstant(1,temp.cols(),axialSliceEnergy);
    avgFissionEnergy.row(iZ) = temporaryMatrix;
  }

  return avgFissionEnergy;

};
//==============================================================================

//==============================================================================
/// Calculate core-average gamma energy deposition term
///
Eigen::MatrixXd HeatTransfer::calcExplicitAxialFuelFissionEnergy()
{
  
  double localVolume,localFlux,totalVolume,localSigF,localOmega,\
    axialSliceVolume = 0,axialSliceEnergy = 0,fuelVol = 0;
 
  Eigen::MatrixXd avgFissionEnergy, temporaryMatrix;
  avgFissionEnergy.setZero(temp.rows(),temp.cols());
 
  totalVolume = M_PI*mesh->R*mesh->R*mesh->Z;
 
  for (int iZ = 0; iZ < temp.rows(); iZ++)
  {
    // Reset slice quantities
    axialSliceEnergy = 0;
    axialSliceVolume = 0;

    for (int iR = 0; iR < temp.cols(); iR++)
    {
      // Get local parameters
      localSigF = mats->oneGroupXS->sigF(iZ,iR);
      localOmega = mats->omega(iZ,iR);
      localFlux = mpqd->ggqd->sFlux(iZ,iR);
      localVolume = mesh->getGeoParams(iR,iZ)[0];

      if (localSigF > 1E-10)
      {
        axialSliceVolume += localVolume;
        axialSliceEnergy += localVolume*(localOmega*localSigF*localFlux); 
      }

    }

    axialSliceEnergy = axialSliceEnergy/axialSliceVolume;
    temporaryMatrix.setConstant(1,temp.cols(),axialSliceEnergy);
    avgFissionEnergy.row(iZ) = temporaryMatrix;
  }

  return avgFissionEnergy;

};
//==============================================================================



//==============================================================================
/// Calculate energy diracs
///
void HeatTransfer::calcDiracs()
{
  double TupwindInterface,Tinterface,theta,phi;
  int lastDiracIndex = dirac.rows()-1;

  if (mats->posVelocity) 
  {
    for (int iR = 0; iR < dirac.cols(); iR++)
    {
      // Handle iZ=0 case
      TupwindInterface = inletTemp(1,iR) - inletTemp(0,iR);
      Tinterface = temp(0,iR) - inletTemp(1,iR);
      theta = calcTheta(TupwindInterface,Tinterface);
      phi = calcPhi(theta,fluxLimiter); 
      dirac(0,iR) = phi*Tinterface;

      // Handle iZ=1 case
      TupwindInterface = temp(0,iR) - inletTemp(1,iR);
      Tinterface = temp(1,iR) - temp(0,iR);
      theta = calcTheta(TupwindInterface,Tinterface);
      phi = calcPhi(theta,fluxLimiter); 
      dirac(1,iR) = phi*Tinterface; 

      // Handle all other cases
      for (int iZ = 2; iZ < dirac.rows()-1; iZ++)
      {

        TupwindInterface = temp(iZ-1,iR) - temp(iZ-2,iR);
        Tinterface = temp(iZ,iR) - temp(iZ-1,iR);
        theta = calcTheta(TupwindInterface,Tinterface);
        phi = calcPhi(theta,fluxLimiter); 
        dirac(iZ,iR) = phi*Tinterface; 

      }

      // Handle iZ = nZ case
      TupwindInterface = temp(lastDiracIndex-1,iR) - temp(lastDiracIndex-2,iR);
      Tinterface = outletTemp(iR) - temp(lastDiracIndex-1,iR);
      theta = calcTheta(TupwindInterface,Tinterface);
      phi = calcPhi(theta,fluxLimiter); 
      dirac(lastDiracIndex,iR) = phi*Tinterface; 

    }
  } else 
  {
    for (int iR = 0; iR < dirac.cols(); iR++)
    {
      // Handle iZ=0 case
      TupwindInterface = temp(1,iR) - temp(0,iR);
      Tinterface = temp(0,iR) - outletTemp(iR);
      theta = calcTheta(TupwindInterface,Tinterface);
      phi = calcPhi(theta,fluxLimiter); 
      dirac(0,iR) = phi*Tinterface; 

      // Handle all other cases
      for (int iZ = 1; iZ < dirac.rows()-2; iZ++)
      {

        TupwindInterface = temp(iZ+1,iR) - temp(iZ,iR);
        Tinterface = temp(iZ,iR) - temp(iZ-1,iR);
        theta = calcTheta(TupwindInterface,Tinterface);
        phi = calcPhi(theta,fluxLimiter); 
        dirac(iZ,iR) = phi*Tinterface; 

      }

      // Handle iZ = nZ-1 case
      TupwindInterface = inletTemp(0,iR) - temp(lastDiracIndex-1,iR);
      Tinterface = temp(lastDiracIndex-1,iR) - temp(lastDiracIndex-2,iR);
      theta = calcTheta(TupwindInterface,Tinterface);
      phi = calcPhi(theta,fluxLimiter); 
      dirac(lastDiracIndex-1,iR) = phi*Tinterface; 

      // Handle iZ = nZ case
      TupwindInterface = inletTemp(1,iR) - inletTemp(0,iR);
      Tinterface = inletTemp(0,iR) - temp(lastDiracIndex-1,iR);
      theta = calcTheta(TupwindInterface,Tinterface);
      phi = calcPhi(theta,fluxLimiter); 
      dirac(lastDiracIndex,iR) = phi*Tinterface; 
    }

  }

};
//==============================================================================

//==============================================================================
/// Calculate energy flux 
///
void HeatTransfer::calcImplicitFluxes()
{

  double tdc; // shorthand for temp*density*specific heat 
  int lastFluxIndex = flux.rows()-1;

  if (mats->posVelocity) {

    for (int iR = 0; iR < flux.cols(); iR++)
    {
      // Handle iZ = 0 case
      tdc = inletVelocity(iR)*inletDensity(iR)*inletcP(iR);
      flux(0,iR) = tdc;

      // Handle all other cases
      for (int iZ = 1; iZ < flux.rows(); iZ++)
      {
        tdc = mats->flowVelocity(iZ-1,iR)*mats->density(iZ-1,iR)\
              *mats->cP(iZ-1,iR);
        flux(iZ,iR) = tdc;

      }

    }

  } else
  {

    for (int iR = 0; iR < flux.cols(); iR++)
    {

      // Handle all other cases
      for (int iZ = 0; iZ < flux.rows()-1; iZ++)
      {
        tdc = mats->flowVelocity(iZ,iR)*mats->density(iZ,iR)\
              *mats->cP(iZ,iR);
        flux(iZ,iR) = tdc;
      }

      // Handle iZ = nZ case
      tdc = inletVelocity(iR)*inletDensity(iR)*inletcP(iR);
      flux(lastFluxIndex,iR) = tdc;

    }
  }

};
//==============================================================================


//==============================================================================
/// Calculate energy flux 
///
void HeatTransfer::calcFluxes()
{

  double tdc; // shorthand for temp*density*specific heat 
  int lastFluxIndex = flux.rows()-1;

  if (mats->posVelocity) {

    for (int iR = 0; iR < flux.cols(); iR++)
    {
      // Handle iZ = 0 case
      tdc = inletVelocity(iR)*inletDensity(iR)*inletcP(iR);
      flux(0,iR) = tdc*inletTemp(1,iR)\
                   + 0.5*abs(tdc)*(1-abs(tdc*mesh->dt/mesh->dzsCorner(0)))\
                   *dirac(0,iR);

      // Handle all other cases
      for (int iZ = 1; iZ < flux.rows(); iZ++)
      {
        tdc = mats->flowVelocity(iZ-1,iR)*mats->density(iZ-1,iR)\
              *mats->cP(iZ-1,iR);
        flux(iZ,iR) = tdc*temp(iZ-1,iR)\
                      + 0.5*abs(tdc)*(1-abs(tdc*mesh->dt/mesh->dzsCorner(iZ-1)))\
                      *dirac(iZ,iR);
      }

    }

  } else
  {

    for (int iR = 0; iR < flux.cols(); iR++)
    {

      // Handle all other cases
      for (int iZ = 0; iZ < flux.rows()-1; iZ++)
      {
        tdc = mats->flowVelocity(iZ,iR)*mats->density(iZ,iR)\
              *mats->cP(iZ,iR);
        flux(iZ,iR) = tdc*temp(iZ,iR)\
                      + 0.5*abs(tdc)*(1-abs(tdc*mesh->dt/mesh->dzsCorner(iZ)))*dirac(iZ,iR);
      }

      // Handle iZ = nZ case
      tdc = inletVelocity(iR)*inletDensity(iR)*inletcP(iR);
      flux(lastFluxIndex,iR) = tdc*inletTemp(0,iR)\
                               + 0.5*abs(tdc)*(1-abs(tdc*mesh->dt/mesh->dzsCorner(lastFluxIndex-1)))\
                               *dirac(lastFluxIndex,iR);

    }
  }

};
//==============================================================================

//==============================================================================
/// Calculate theta for use in flux limiting framework 
///
/// @param [in] TupwindInterface change in temp at upwind interface
/// @param [in] Tinterface change in temp at interface
double HeatTransfer::calcTheta(double TupwindInterface, double Tinterface)
{

  double theta;
  if (abs(Tinterface) < 1E-10)
  {
    theta = 1;
  } else 
  {
    theta = TupwindInterface/Tinterface;
  }

  return theta;

};
//==============================================================================


//==============================================================================
/// Calculate phi for use in flux limiting framework 
///
/// @param [in] theta ratio of change in temperature in upwind and current cell
/// @param [in] fluxLimiter string indicating how to calculate phi
/// @param [out] phi flux limiting parameter
double HeatTransfer::calcPhi(double theta, string fluxLimiter)
{

  double phi; 
  Eigen::Vector2d fluxLimiterArg1,fluxLimiterArg2;
  Eigen::Vector3d fluxLimiterArg3;

  if (fluxLimiter == "superbee")
  {

    fluxLimiterArg1 << 1,2*theta; 
    fluxLimiterArg2 << 2,theta; 
    fluxLimiterArg3 << 0,\
      fluxLimiterArg1.minCoeff(),\
      fluxLimiterArg2.minCoeff();
    phi = fluxLimiterArg3.maxCoeff();

  } else if (fluxLimiter == "upwind") 
  {
    phi = 0.0;
  } else if (fluxLimiter == "lax-wendroff")
  {
    phi = 1.0;
  } else if (fluxLimiter == "beam-warming")
  {
    phi = theta;
  }

  return phi;

};
//==============================================================================

//==============================================================================
/// Assign boundary indices depending on direction of flow velocity
///
void HeatTransfer::assignBoundaryIndices()
{
  if (mats->posVelocity)
  { 
    coreInletIndex = 0;
    coreOutletIndex = mesh->nZ-1;
  } else 
  {
    coreInletIndex = mesh->nZ-1;
    coreOutletIndex = 0;
  }

};
//==============================================================================

//==============================================================================
/// Update temperatures and physical parameters at inlet and outlet     
///
void HeatTransfer::updateBoundaryConditions()
{

  // Update variables that can set with array splicing
  inletTemp.setConstant(inletTemp.rows(),inletTemp.cols(),inletT);
  outletTemp = temp.row(coreOutletIndex);   
  inletVelocity = mats->flowVelocity.row(coreInletIndex);

  // Update variables that require looping to access   
  for (int iR = 0; iR < mesh->nR; iR++)
  {
    inletDensity(iR) = mats->density(coreInletIndex,iR);
    inletcP(iR) = mats->cP(coreInletIndex,iR);
  }

};
//==============================================================================

//==============================================================================
/// Map 2D coordinates to index of temperature in the 1D solution vector
///
/// @param [in] iZ axial location
/// @param [in] iR radial location
/// @param [out] index the index for temperature in the 1D solution vector
int HeatTransfer::getIndex(int iZ, int iR)
{

  int index,nR = mesh->drsCorner.size();

  index = indexOffset + iR + nR*iZ;

  return index;

};
//==============================================================================

//==============================================================================
/// Map 2D coordinates to index of temperature in the 1D solution vector
///
void HeatTransfer::getTemp()
{

  if (not mats->uniformTemperature)
    temp = returnCurrentTemp();
  else
    temp.setConstant(mats->uniformTempValue);


};
//==============================================================================

//==============================================================================
/// Map 2D coordinates to index of temperature in the 1D solution vector
///
Eigen::MatrixXd HeatTransfer::returnCurrentTemp()
{
  Eigen::MatrixXd currentTemp; 
  PetscScalar value[1]; 
  PetscInt index[1]; 
  VecScatter     ctx;
  Vec temp_x_p_seq;

  currentTemp.setZero(mesh->nZ,mesh->nR);

  if (mesh->petsc)
  {

    // Gather values of x_p on all procs
    VecScatterCreateToAll(mpqd->x_p,&ctx,&temp_x_p_seq);
    VecScatterBegin(ctx,mpqd->x_p,temp_x_p_seq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,mpqd->x_p,temp_x_p_seq,INSERT_VALUES,SCATTER_FORWARD);

    // loop over spatial mesh
    for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
    {
      for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
      {

        // Read fluxes into flux vector
        index[0] = getIndex(iZ,iR);
        VecGetValues(temp_x_p_seq,1,index,value);
        currentTemp(iZ,iR) = value[0];   

      }
    } 
        
    /* Destroy scatter context */
    VecScatterDestroy(&ctx);
    VecDestroy(&temp_x_p_seq);

  }
  else
  {
    for (int iR = 0; iR < mesh-> drsCorner.size(); iR++)
    {
      for (int iZ = 0; iZ < mesh-> dzsCorner.size(); iZ++)
      { 

        currentTemp(iZ,iR) = mpqd->x(getIndex(iZ,iR));   

      }
    }
  }

  return currentTemp; 
};
//==============================================================================

//==============================================================================
/// Map values from 2D matrix into 1D solution vector
///
int HeatTransfer::setTemp()
{

  PetscErrorCode ierr;
  PetscScalar value;
  PetscInt index;

  if (mesh->petsc)
  {
    for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
    {
      for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
      { 

        value = temp(iZ,iR);
        index = getIndex(iZ,iR);
        ierr = VecSetValue(mpqd->xPast_p,index,value,ADD_VALUES);CHKERRQ(ierr); 

      }
    }
    
    ierr = VecAssemblyBegin(mpqd->xPast_p);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(mpqd->xPast_p);CHKERRQ(ierr);
  }
  else
  {
    for (int iR = 0; iR < mesh-> drsCorner.size(); iR++)
    {
      for (int iZ = 0; iZ < mesh-> dzsCorner.size(); iZ++)
      { 

        mpqd->xPast(getIndex(iZ,iR)) = temp(iZ,iR);   

      }
    }
  }

  return ierr;

};
//==============================================================================

/* PETSc functions */

// Steady state

//==============================================================================
/// Build steady state linear system to solve heat equation
///
int HeatTransfer::buildSteadyStateLinearSystem_p()
{

  int myIndex,sIndex,nIndex,wIndex,eIndex,upwindIndex,iEq = indexOffset;
  int iEqTemp = 0;
  int nR = temp.cols()-1;
  int nZ = temp.rows()-1;
  double harmonicAvg,coeff,keff,neutronFlux,cCoeff;
  Eigen::MatrixXd volAvgGammaDep;
  vector<double> gParams;
  PetscErrorCode ierr;
  PetscScalar value;

  updateBoundaryConditions();
  calcImplicitFluxes();

  // Calculate core-average gamma deposition term
  if (modIrradiation == axial)
    volAvgGammaDep = calcExplicitAxialFissionEnergy();
  else if (modIrradiation == fuel)
    volAvgGammaDep = calcExplicitAxialFuelFissionEnergy();
  else
    volAvgGammaDep = calcExplicitFissionEnergy();

  //#pragma omp parallel for private(myIndex,sIndex,nIndex,wIndex,eIndex,\
  upwindIndex,gParams,cCoeff,coeff,keff,neutronFlux,harmonicAvg,iEq,iEqTemp)
    for (int iZ = 0; iZ < temp.rows(); iZ++)
    {
      for (int iR = 0; iR < temp.cols(); iR++)
      {

        iEq = getIndex(iZ,iR);
        iEqTemp = iEq - indexOffset;

        // Reset center coefficient
        cCoeff = 0;

        // Get cell indices
        myIndex = getIndex(iZ,iR);     
        sIndex = getIndex(iZ+1,iR);     
        nIndex = getIndex(iZ-1,iR);     
        wIndex = getIndex(iZ,iR-1);     
        eIndex = getIndex(iZ,iR+1);     

        gParams = mesh->getGeoParams(iR,iZ);

        // East face
        if (iR == nR)
        {
          coeff = (-gParams[iEF]*mats->k(iZ,iR)\
              /mesh->drsCorner(iR))/gParams[iVol];
          cCoeff -= coeff;
          value = -coeff*wallT;
          ierr = VecSetValue(mpqd->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
          //mpqd->b(iEq) -= coeff*wallT; 
        } else
        {
          harmonicAvg = pow(mesh->drsCorner(iR)/mats->k(iZ,iR)\
              + mesh->drsCorner(iR+1)/mats->k(iZ,iR+1),-1.0);
          coeff = -2.0*gParams[iEF]*harmonicAvg/gParams[iVol];
          //Atemp(iEqTemp,eIndex) = coeff;
          ierr = MatSetValue(mpqd->A_p,iEq,eIndex,coeff,ADD_VALUES);CHKERRQ(ierr); 
          cCoeff -= coeff;
        }

        // West face
        if (iR != 0)
        {
          harmonicAvg = pow(mesh->drsCorner(iR-1)/mats->k(iZ,iR-1)\
              + mesh->drsCorner(iR)/mats->k(iZ,iR),-1.0);
          coeff = 2.0*gParams[iWF]*harmonicAvg/gParams[iVol];
          ierr = MatSetValue(mpqd->A_p,iEq,wIndex,-coeff,ADD_VALUES);CHKERRQ(ierr); 
          //Atemp(iEqTemp,wIndex) = -coeff;
          cCoeff += coeff;
        } 

        // North face
        if (iZ == 0 and mats->posVelocity)
        {
          coeff = (gParams[iNF]*mats->k(iZ,iR)\
              /mesh->dzsCorner(iZ))/gParams[iVol];
          cCoeff += coeff;
          value = coeff*inletTemp(1,iR);
          ierr = VecSetValue(mpqd->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
          //mpqd->b(iEq) += coeff*inletTemp(1,iR);             
        } else if (iZ != 0)
        {
          harmonicAvg = pow(mesh->dzsCorner(iZ-1)/mats->k(iZ-1,iR)\
              + mesh->dzsCorner(iZ)/mats->k(iZ,iR),-1.0);
          coeff = 2.0*gParams[iNF]*harmonicAvg/gParams[iVol];
          ierr = MatSetValue(mpqd->A_p,iEq,nIndex,-coeff,ADD_VALUES);CHKERRQ(ierr); 
          //Atemp(iEqTemp,nIndex) = -coeff;
          cCoeff += coeff;
        }

        // South face
        if (iZ == nZ and !(mats->posVelocity))
        {
          coeff = -(gParams[iSF]*mats->k(iZ,iR)\
              /mesh->dzsCorner(iZ))/gParams[iVol];
          cCoeff -= coeff;
          value = coeff*inletTemp(0,iR);
          ierr = VecSetValue(mpqd->b_p,iEq,-value,ADD_VALUES);CHKERRQ(ierr); 
          //mpqd->b(iEq) -= coeff*inletTemp(0,iR);             
        } else if (iZ != nZ)
        {
          harmonicAvg = pow(mesh->dzsCorner(iZ+1)/mats->k(iZ+1,iR)\
              + mesh->dzsCorner(iZ)/mats->k(iZ,iR),-1.0);
          coeff = -2.0*gParams[iSF]*harmonicAvg/gParams[iVol];
          ierr = MatSetValue(mpqd->A_p,iEq,sIndex,coeff,ADD_VALUES);CHKERRQ(ierr); 
          //Atemp(iEqTemp,sIndex) = coeff;
          cCoeff -= coeff;
        }

        // Insert cell center coefficient
        //Atemp.insert(iEqTemp,myIndex) = cCoeff;
        //Atemp(iEqTemp,myIndex) = cCoeff;
        ierr = MatSetValue(mpqd->A_p,iEq,myIndex,cCoeff,ADD_VALUES);CHKERRQ(ierr); 

        // Flux source term 
        coeff = mats->omega(iZ,iR)*mats->oneGroupXS->sigF(iZ,iR);
        keff = mats->oneGroupXS->keff; 
        neutronFlux = mpqd->ggqd->sFlux(iZ,iR);

        //mpqd->b(iEq) += coeff*neutronFlux/keff; 
        value = coeff*neutronFlux;
        ierr = VecSetValue(mpqd->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
        //mpqd->b(iEq) += coeff*neutronFlux; 

        // Gamma source term 
        value = mats->gamma(iZ,iR) * volAvgGammaDep(iZ,iR);
        ierr = VecSetValue(mpqd->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
        //mpqd->b(iEq) += mats->gamma(iZ,iR) * volAvgGammaDep(iZ,iR);

        if (mats->posVelocity) 
        {
          upwindIndex = getIndex(iZ-1,iR);     

          // Advection terms

          // Upwind cell
          if (iZ == 0) // boundary case
          {
            value = flux(iZ,iR)*inletTemp(1,iR)/mesh->dzsCorner(iZ);
            ierr = VecSetValue(mpqd->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
            //mpqd->b(iEq) += flux(iZ,iR)*inletTemp(1,iR)/mesh->dzsCorner(iZ);
          }
          else
          {
            value = -flux(iZ,iR)/mesh->dzsCorner(iZ);
            ierr = MatSetValue(mpqd->A_p,iEq,upwindIndex,value,ADD_VALUES);CHKERRQ(ierr); 
            //Atemp.coeffRef(iEqTemp,upwindIndex) += -flux(iZ,iR)/mesh->dzsCorner(iZ);
          }

          // Primary cell
          value = flux(iZ+1,iR)/mesh->dzsCorner(iZ);
          ierr = MatSetValue(mpqd->A_p,iEq,myIndex,value,ADD_VALUES);CHKERRQ(ierr); 
          //Atemp.coeffRef(iEqTemp,myIndex) += flux(iZ+1,iR)/mesh->dzsCorner(iZ);
        }
        else
        {
          upwindIndex = getIndex(iZ+1,iR);     

          // Advection terms

          // Upwind cell
          if (iZ == temp.rows()) // boundary case
          {
            value = flux(iZ+1,iR)*inletTemp(0,iR)/mesh->dzsCorner(iZ);
            ierr = VecSetValue(mpqd->b_p,iEq,-value,ADD_VALUES);CHKERRQ(ierr); 
            //mpqd->b(iEq) -= flux(iZ+1,iR)*inletTemp(0,iR)/mesh->dzsCorner(iZ);
          }
          else
          {
            value = flux(iZ+1,iR)/mesh->dzsCorner(iZ);
            ierr = MatSetValue(mpqd->A_p,iEq,upwindIndex,value,ADD_VALUES);CHKERRQ(ierr); 
            //Atemp.coeffRef(iEqTemp,upwindIndex) = flux(iZ+1,iR)/mesh->dzsCorner(iZ);
          }

          // Primary cell
          value = flux(iZ,iR)/mesh->dzsCorner(iZ);
          ierr = MatSetValue(mpqd->A_p,iEq,myIndex,-value,ADD_VALUES);CHKERRQ(ierr); 
          //Atemp.coeffRef(iEqTemp,myIndex) -= flux(iZ,iR)/mesh->dzsCorner(iZ);
        }
      }
    }

  //mpqd->A.middleRows(indexOffset,nUnknowns) = Atemp.sparseView(); 

  return ierr;

};
//==============================================================================

// TRANSIENT

//==============================================================================
/// Build linear system to solve heat equation
///
int HeatTransfer::buildLinearSystem_p()
{
  
  int myIndex,sIndex,nIndex,wIndex,eIndex,iEq = indexOffset;
  int iEqTemp = 0;
  int nR = temp.cols()-1;
  int nZ = temp.rows()-1;
  double harmonicAvg,coeff,cCoeff;
  Eigen::MatrixXd volAvgGammaDep;
  vector<double> gParams;
  PetscErrorCode ierr;
  PetscScalar value;

  updateBoundaryConditions();
  calcDiracs();
  calcFluxes();

  // Calculate core-average gamma deposition term
  if (modIrradiation == axial)
    volAvgGammaDep = calcExplicitAxialFissionEnergy();
  else if (modIrradiation == fuel)
    volAvgGammaDep = calcExplicitAxialFuelFissionEnergy();
  else
    volAvgGammaDep = calcExplicitFissionEnergy();

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Atemp;
  
  //#pragma omp parallel for private(myIndex,sIndex,nIndex,wIndex,eIndex,\
    gParams,cCoeff,coeff,harmonicAvg,iEq,iEqTemp)
  for (int iZ = 0; iZ < temp.rows(); iZ++)
  {

    for (int iR = 0; iR < temp.cols(); iR++)
    {

      iEq = getIndex(iZ,iR);
      iEqTemp = iEq - indexOffset;

      // Reset center coefficient
      cCoeff = 0;

      // Get cell indices
      myIndex = getIndex(iZ,iR);     
      sIndex = getIndex(iZ+1,iR);     
      nIndex = getIndex(iZ-1,iR);     
      wIndex = getIndex(iZ,iR-1);     
      eIndex = getIndex(iZ,iR+1);     

      gParams = mesh->getGeoParams(iR,iZ);
      
      cCoeff = mats->density(iZ,iR)*mats->cP(iZ,iR);
      
      // East face
      if (iR == nR)
      {
        coeff = (-gParams[iEF]*mats->k(iZ,iR)\
        /mesh->drsCorner(iR))/gParams[iVol];
        value = -mesh->dt*coeff*wallT;
        ierr = VecSetValue(mpqd->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
        cCoeff -= mesh->dt*coeff;
        //mpqd->b(iEq) -= mesh->dt*coeff*wallT; 
      } else
      {
        harmonicAvg = pow(mesh->drsCorner(iR)/mats->k(iZ,iR)\
          + mesh->drsCorner(iR+1)/mats->k(iZ,iR+1),-1.0);
        coeff = -2.0*gParams[iEF]*harmonicAvg/gParams[iVol];
        value = mesh->dt*coeff;
        ierr = MatSetValue(mpqd->A_p,iEq,eIndex,value,ADD_VALUES);CHKERRQ(ierr); 
        cCoeff -= mesh->dt*coeff;
        //Atemp(iEqTemp,eIndex) = mesh->dt*coeff;
      }

      // West face
      if (iR != 0)
      {
        harmonicAvg = pow(mesh->drsCorner(iR-1)/mats->k(iZ,iR-1)\
          + mesh->drsCorner(iR)/mats->k(iZ,iR),-1.0);
        coeff = 2.0*gParams[iWF]*harmonicAvg/gParams[iVol];
        value = -mesh->dt*coeff;
        ierr = MatSetValue(mpqd->A_p,iEq,wIndex,value,ADD_VALUES);CHKERRQ(ierr); 
        cCoeff += mesh->dt*coeff;
        //Atemp(iEqTemp,wIndex) = -mesh->dt*coeff;
      } 

      // North face
      if (iZ == 0 and mats->posVelocity)
      {
        coeff = (gParams[iNF]*mats->k(iZ,iR)\
        /mesh->dzsCorner(iZ))/gParams[iVol];
        cCoeff += mesh->dt*coeff;
        value = mesh->dt*coeff*inletTemp(1,iR);
        ierr = VecSetValue(mpqd->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
        //mpqd->b(iEq) += mesh->dt*coeff*inletTemp(1,iR);             
      } else if (iZ != 0)
      {
        harmonicAvg = pow(mesh->dzsCorner(iZ-1)/mats->k(iZ-1,iR)\
          + mesh->dzsCorner(iZ)/mats->k(iZ,iR),-1.0);
        coeff = 2.0*gParams[iNF]*harmonicAvg/gParams[iVol];
        cCoeff += mesh->dt*coeff;
        value = -mesh->dt*coeff;
        ierr = MatSetValue(mpqd->A_p,iEq,nIndex,value,ADD_VALUES);CHKERRQ(ierr); 
        //Atemp(iEqTemp,nIndex) = -mesh->dt*coeff;
      }

      // South face
      if (iZ == nZ and !(mats->posVelocity))
      {
        coeff = -(gParams[iSF]*mats->k(iZ,iR)\
        /mesh->dzsCorner(iZ))/gParams[iVol];
        cCoeff -= mesh->dt*coeff;
        value = -mesh->dt*coeff*inletTemp(0,iR);
        ierr = VecSetValue(mpqd->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
        //mpqd->b(iEq) -= mesh->dt*coeff*inletTemp(0,iR);             
      } else if (iZ != nZ)
      {
        harmonicAvg = pow(mesh->dzsCorner(iZ+1)/mats->k(iZ+1,iR)\
          + mesh->dzsCorner(iZ)/mats->k(iZ,iR),-1.0);
        coeff = -2.0*gParams[iSF]*harmonicAvg/gParams[iVol];
        value = mesh->dt*coeff;
        ierr = MatSetValue(mpqd->A_p,iEq,sIndex,value,ADD_VALUES);CHKERRQ(ierr); 
        cCoeff -= mesh->dt*coeff;
        //Atemp(iEqTemp,sIndex) = mesh->dt*coeff;
      }

      // Insert cell center coefficient
      ierr = MatSetValue(mpqd->A_p,iEq,myIndex,cCoeff,ADD_VALUES);CHKERRQ(ierr); 
      //Atemp(iEqTemp,myIndex) = cCoeff;

      // Time term
      value = mats->density(iZ,iR)*mats->cP(iZ,iR)*temp(iZ,iR); 
      ierr = VecSetValue(mpqd->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
      //mpqd->b(iEq) += mats->density(iZ,iR)*mats->cP(iZ,iR)*temp(iZ,iR); 

      // Flux source term 
      coeff = -mesh->dt*mats->omega(iZ,iR)*mats->oneGroupXS->sigF(iZ,iR);
      mpqd->fluxSource(iZ,iR,iEq,coeff,&Atemp);
      
      // Gamma source term 
      coeff = mesh->dt;
      value = coeff * mats->gamma(iZ,iR) * volAvgGammaDep(iZ,iR);
      ierr = VecSetValue(mpqd->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
      //mpqd->b(iEq) += coeff * mats->gamma(iZ,iR) * volAvgGammaDep(iZ,iR);
      //gammaSource(iZ,iR,iEqTemp,coeff);

      // Advection term
      value = (mesh->dt/mesh->dzsCorner(iZ))*(flux(iZ,iR)-flux(iZ+1,iR));
      ierr = VecSetValue(mpqd->b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
      //mpqd->b(iEq) += (mesh->dt/mesh->dzsCorner(iZ))*(flux(iZ,iR)-flux(iZ+1,iR));

      // Iterate equation count
      //iEq = iEq + 1;
      //iEqTemp = iEqTemp + 1;
  
    }
  }

  return ierr;

};
//==============================================================================

//==============================================================================
/// Check for optional input parameters of relevance to this object
void HeatTransfer::checkOptionalParams()
{

  // Check for optional inputs 
  if ((*input)["parameters"]["wallTemp"]){
    wallT=(*input)["parameters"]["wallTemp"].as<double>();
  } 
  if ((*input)["parameters"]["inletTemp"]){
    inletT=(*input)["parameters"]["inletTemp"].as<double>();
  } 
  if ((*input)["parameters"]["flux limiter"]){
    fluxLimiter=(*input)["parameters"]["flux limiter"].as<string>();
  } 
  if ((*input)["parameters"]["modIrradiation"])
  {
    modIrradiation=(*input)["parameters"]["modIrradiation"].as<string>();
  }

}
//==============================================================================
