// File: SingleGroupDNP.cpp     
// Purpose: Contains information and builds linear system for a precursor group
// Date: April 9, 2020

#include "SingleGroupDNP.h"
#include "MultiGroupDNP.h"
#include "MultiPhysicsCoupledQD.h"

using namespace std;

//==============================================================================
/// SingleGroupDNP class object constructor
///
/// @param [in] blankType blank for this material
SingleGroupDNP::SingleGroupDNP(Materials * myMats,\
  Mesh * myMesh,\
  MultiGroupDNP * myMGDNP,\
  double myBeta,\
  double myLambda,\
  int myCoreIndexOffset,\
  int myRecircIndexOffset)
{

  // Assign inputs to their member variables
  mats = myMats;
  mesh = myMesh;
  mgdnp = myMGDNP;
  beta = myBeta;
  lambda = myLambda;
  coreIndexOffset = myCoreIndexOffset;
  recircIndexOffset = myRecircIndexOffset;

  // Initialize size of matrices and vectors
  dnpConc.setOnes(mesh->nZ,mesh->nR);
  flux.setZero(mesh->nZ+1,mesh->nR); 
  dirac.setZero(mesh->nZ+1,mesh->nR);
  recircConc.setOnes(mesh->nZrecirc,mesh->nR); 
  recircFlux.setZero(mesh->nZrecirc+1,mesh->nR); 
  recircDirac.setZero(mesh->nZrecirc+1,mesh->nR);
  inletConc.setZero(2,mesh->nR);
  recircInletConc.setZero(2,mesh->nR);
  inletVelocity.setZero(mesh->nR);
  recircInletVelocity.setZero(mesh->nR);
  outletConc.setZero(mesh->nR);
  recircOutletConc.setZero(mesh->nR);

  // assign boundary conditions depending on direction of flow
  assignBoundaryIndices();
 
};
//==============================================================================

//==============================================================================
/// Build linear system for this precursor group
///
void SingleGroupDNP::buildLinearSystem(Eigen::SparseMatrix<double> * myA,\
  Eigen::VectorXd * myb,\
  Eigen::MatrixXd myDNPConc,\
  Eigen::MatrixXd myDNPFlux,\
  Eigen::MatrixXd mySigF,\
  rowvec dzs,\
  int myIndexOffset,\
  bool fluxSource)
{

  int myIndex,iEq = myIndexOffset;
  double coeff;
  
  for (int iR = 0; iR < myDNPConc.cols(); iR++)
  {
    for (int iZ = 0; iZ < myDNPConc.rows(); iZ++)
    {
      myIndex = getIndex(iZ,iR,myIndexOffset);     

      myA->coeffRef(iEq,myIndex) = 1 + mesh->dt*lambda; 
      
      // Time term
      (*myb)(iEq) = myDNPConc(iZ,iR);

      // Flux source term 
      coeff = mesh->dt*beta*mySigF(iZ,iR);
      
      if (fluxSource)
        mgdnp->mpqd->fluxSource(iZ,iR,iEq,coeff);

      // Advection term
      (*myb)(iEq) += (mesh->dt/dzs(iZ))*(myDNPFlux(iZ,iR)-myDNPFlux(iZ+1,iR));

      // Iterate equation count
      iEq = iEq + 1;
  
    }
  }

};
//==============================================================================


//==============================================================================
/// Calculate diracs to model advection of precursors
///
Eigen::MatrixXd SingleGroupDNP::calcDiracs(Eigen::MatrixXd dnpConc,\
  Eigen::MatrixXd inletConc,\
  Eigen::VectorXd outletConc)
{

  Eigen::MatrixXd myDirac(dnpConc.rows()+1,dnpConc.cols());
  double CupwindInterface,Cinterface,theta,phi;
  int lastDiracIndex = myDirac.rows()-1;

  if (mats->posVelocity) 
  {
    for (int iR = 0; iR < myDirac.cols(); iR++)
    {
      // Handle iZ=0 case
      CupwindInterface = inletConc(1,iR) - inletConc(0,iR);
      Cinterface = dnpConc(0,iR) - inletConc(1,iR);
      theta = calcTheta(CupwindInterface,Cinterface);
      phi = calcPhi(theta,fluxLimiter); 
      myDirac(0,iR) = phi*Cinterface;

      // Handle iZ=1 case
      CupwindInterface = dnpConc(0,iR) - inletConc(1,iR);
      Cinterface = dnpConc(1,iR) - dnpConc(0,iR);
      theta = calcTheta(CupwindInterface,Cinterface);
      phi = calcPhi(theta,fluxLimiter); 
      myDirac(1,iR) = phi*Cinterface; 
      
      // Handle all other cases
      for (int iZ = 2; iZ < myDirac.rows()-1; iZ++)
      {

        CupwindInterface = dnpConc(iZ-1,iR) - dnpConc(iZ-2,iR);
        Cinterface = dnpConc(iZ,iR) - dnpConc(iZ-1,iR);
        theta = calcTheta(CupwindInterface,Cinterface);
        phi = calcPhi(theta,fluxLimiter); 
        myDirac(iZ,iR) = phi*Cinterface; 
        
      }

      // Handle iZ = nZ case
      CupwindInterface = dnpConc(lastDiracIndex-1,iR)\
        - dnpConc(lastDiracIndex-2,iR);
      Cinterface = outletConc(iR) - dnpConc(lastDiracIndex-1,iR);
      theta = calcTheta(CupwindInterface,Cinterface);
      phi = calcPhi(theta,fluxLimiter); 
      myDirac(lastDiracIndex,iR) = phi*Cinterface;
      
    }
  } else 
  {
    for (int iR = 0; iR < myDirac.cols(); iR++)
    {
      // Handle iZ=0 case
      CupwindInterface = dnpConc(1,iR) - dnpConc(0,iR);
      Cinterface = dnpConc(0,iR) - outletConc(iR);
      theta = calcTheta(CupwindInterface,Cinterface);
      phi = calcPhi(theta,fluxLimiter); 
      myDirac(0,iR) = phi*Cinterface; 
      
      // Handle all other cases
      for (int iZ = 1; iZ < myDirac.rows()-2; iZ++)
      {

        CupwindInterface = dnpConc(iZ+1,iR) - dnpConc(iZ,iR);
        Cinterface = dnpConc(iZ,iR) - dnpConc(iZ-1,iR);
        theta = calcTheta(CupwindInterface,Cinterface);
        phi = calcPhi(theta,fluxLimiter); 
        myDirac(iZ,iR) = phi*Cinterface; 
        
      }

      // Handle iZ = nZ-1 case
      CupwindInterface = inletConc(0,iR) - dnpConc(lastDiracIndex-1,iR);
      Cinterface = dnpConc(lastDiracIndex-1,iR) - dnpConc(lastDiracIndex-2,iR);
      theta = calcTheta(CupwindInterface,Cinterface);
      phi = calcPhi(theta,fluxLimiter); 
      myDirac(lastDiracIndex-1,iR) = phi*Cinterface; 

      // Handle iZ = nZ case
      CupwindInterface = inletConc(1,iR) - inletConc(0,iR);
      Cinterface = inletConc(0,iR) - dnpConc(lastDiracIndex-1,iR);
      theta = calcTheta(CupwindInterface,Cinterface);
      phi = calcPhi(theta,fluxLimiter); 
      myDirac(lastDiracIndex,iR) = phi*Cinterface; 
    }

  }
  
  cout << "calculated diracs" << endl;
  cout << myDirac << endl;

  return myDirac;

};
//==============================================================================

//==============================================================================
/// Calculate fluxes to model advection of precursors
///
Eigen::MatrixXd SingleGroupDNP::calcFluxes(Eigen::MatrixXd myDNPConc,\
  Eigen::MatrixXd myFlowVelocity,\
  Eigen::MatrixXd myDirac,\
  Eigen::MatrixXd myInletConc,\
  Eigen::VectorXd myInletVelocity,\
  rowvec dzs)
{

  Eigen::MatrixXd myFlux;
  double vel; // shorthand for velocity
  int lastFluxIndex = myDirac.rows()-1;

  myFlux.setZero(myDirac.rows(),myDirac.cols());
 
  if (mats->posVelocity) {

    for (int iR = 0; iR < myFlux.cols(); iR++)
    {
      // Handle iZ = 0 case
      vel = myInletVelocity(iR);
      myFlux(0,iR) = vel*myInletConc(1,iR)\
        + 0.5*abs(vel)*(1-abs(vel*mesh->dt/dzs(0)))\
        *myDirac(0,iR);

      // Handle all other cases
      for (int iZ = 1; iZ < myFlux.rows(); iZ++)
      {
        vel = myFlowVelocity(iZ-1,iR);
        myFlux(iZ,iR) = vel*myDNPConc(iZ-1,iR)\
          + 0.5*abs(vel)*(1-abs(vel*mesh->dt/dzs(iZ-1)))\
          *myDirac(iZ,iR);
      }

    }
    
  } else
  {

    for (int iR = 0; iR < myFlux.cols(); iR++)
    {

      // Handle all other cases
      for (int iZ = 0; iZ < myFlux.rows()-1; iZ++)
      {
        vel = myFlowVelocity(iZ,iR);
        myFlux(iZ,iR) = vel*myDNPConc(iZ,iR)\
          + 0.5*abs(vel)*(1-abs(vel*mesh->dt/dzs(iZ)))*myDirac(iZ,iR);
      }
      
      // Handle iZ = nZ case
      vel = myInletVelocity(iR);
      myFlux(lastFluxIndex,iR) = vel*myInletConc(0,iR)\
        + 0.5*abs(vel)*(1-abs(vel*mesh->dt/dzs(lastFluxIndex-1)))\
        *myDirac(lastFluxIndex,iR);

    }
  }
  
  cout << "calculated fluxes" << endl;
  cout << myFlux << endl;

  return myFlux;

};
//==============================================================================

//==============================================================================
/// Map 2D coordinates to index of delayed neutron precursor concentration in 
/// the 1D solution
///
/// @param [in] iZ axial location
/// @param [in] iR radial location
/// @param [out] index the index for dnp concentration in the 1D solution vector
int SingleGroupDNP::getIndex(int iZ, int iR,int indexOffset)
{

  int index,nR = mesh->drsCorner.size();

  index = indexOffset + iR + nR*iZ;

  return index;
  
};
//==============================================================================

//==============================================================================
/// Map solution in 1D vector to 2D solution
///
void SingleGroupDNP::getCoreConc()
{

  for (int iR = 0; iR < mesh-> drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh-> dzsCorner.size(); iZ++)
    { 
    
      dnpConc(iZ,iR) = mgdnp->mpqd->x(getIndex(iZ,iR,coreIndexOffset));   
    
    }
  }
  
};
//==============================================================================

//==============================================================================
/// Map solution in 1D vector to 2D solution
///
void SingleGroupDNP::getRecircConc()
{

  for (int iR = 0; iR < mesh-> drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh-> nZrecirc; iZ++)
    { 
    
      recircConc(iZ,iR) = mgdnp->recircx(getIndex(iZ,iR,recircIndexOffset));   
    
    }
  }
  
};
//==============================================================================


//==============================================================================
/// Calculate phi factor in flux limiting scheme
///
/// @param [in] theta ratio of change in temperature in upwind and current cell
/// @param [in] fluxLimiter string indicating how to calculate phi
/// @param [out] phi flux limiting parameter
double SingleGroupDNP::calcPhi(double theta,string fluxLimiter)
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
/// Calculate theta factor in flux limiting scheme
///
/// @param [in] DNPupwindInterface change in temp at upwind interface
/// @param [in] DNPinterface change in temp at interface
/// @param [out] theta ratio of deltas at interfaces
double SingleGroupDNP::calcTheta(double DNPupwindInterface,double DNPinterface)
{

  double theta;
  if (abs(DNPinterface) < 1E-10)
  {
    theta = 1;
  } else 
  {
    theta = DNPupwindInterface/DNPinterface;
  }

  return theta;

};
//==============================================================================

//==============================================================================
/// Calculate theta factor in flux limiting scheme
///
void SingleGroupDNP::calcRecircDNPFluxes()
{

  Eigen::MatrixXd recircDirac,recircFlux;

  recircDirac = calcDiracs(recircConc,\
    recircInletConc,\
    recircOutletConc);

  recircFlux = calcFluxes(recircConc,\
    mats->recircFlowVelocity,\
    recircDirac,\
    recircInletConc,\
    recircInletVelocity,\
    mesh->dzsCornerRecirc);

};
//==============================================================================

//==============================================================================
/// Calculate theta factor in flux limiting scheme
///
void SingleGroupDNP::calcCoreDNPFluxes()
{

  Eigen::MatrixXd coreDirac,coreFlux;

  coreDirac = calcDiracs(dnpConc,\
    inletConc,\
    outletConc);

  coreFlux = calcFluxes(dnpConc,\
    mats->flowVelocity,\
    coreDirac,\
    inletConc,\
    inletVelocity,\
    mesh->dzsCorner);

};
//==============================================================================

//==============================================================================
/// Calculate theta factor in flux limiting scheme
///
void SingleGroupDNP::buildRecircLinearSystem()
{

  Eigen::MatrixXd recircDirac,recircFlux,dumbySigF;

  // Assume there are no source terms in the recirculation loop 
  dumbySigF.setZero(mesh->nZrecirc,mesh->nR);

  updateBoundaryConditions();

  recircDirac = calcDiracs(recircConc,\
    recircInletConc,\
    recircOutletConc);

  recircFlux = calcFluxes(recircConc,\
    mats->recircFlowVelocity,\
    recircDirac,\
    recircInletConc,\
    recircInletVelocity,\
    mesh->dzsCornerRecirc);

  buildLinearSystem(&(mgdnp->recircA),\
    &(mgdnp->recircb),\
    recircConc,\
    recircFlux,\
    dumbySigF,\
    mesh->dzsCornerRecirc,\
    recircIndexOffset,\
    false);

};
//==============================================================================

//==============================================================================
/// Calculate theta factor in flux limiting scheme
///
void SingleGroupDNP::buildCoreLinearSystem()
{

  Eigen::MatrixXd coreDirac,coreFlux;
  
  updateBoundaryConditions();

  coreDirac = calcDiracs(dnpConc,\
    inletConc,\
    outletConc);

  coreFlux = calcFluxes(dnpConc,\
    mats->flowVelocity,\
    coreDirac,\
    inletConc,\
    inletVelocity,\
    mesh->dzsCorner);

  buildLinearSystem(&(mgdnp->mpqd->A),\
    &(mgdnp->mpqd->b),\
    dnpConc,\
    coreFlux,\
    mats->oneGroupXS->sigF,\
    mesh->dzsCorner,\
    coreIndexOffset);
};
//==============================================================================

//==============================================================================
/// Assign boundary indices depending on direction of flow velocity
///
void SingleGroupDNP::assignBoundaryIndices()
{

  if (mats->posVelocity)
  {
    // Assign core indices
    coreInletIndex = mesh->nZrecirc-1;
    coreOutletIndex = 0;

    // Assign recirculation indices
    recircInletIndex = mesh->nZ-1;
    recircOutletIndex = 0;

  } else
  {
    // Assign core indices
    coreInletIndex = 0;
    coreOutletIndex = mesh->nZrecirc-1;
    
    // Assign recirculation indices
    recircInletIndex = 0;
    recircOutletIndex = mesh->nZ-1;
  }
 
};
//==============================================================================

//==============================================================================
/// Update boundary conditions 
///
void SingleGroupDNP::updateBoundaryConditions()
{
  cout << "inletConc:"<<inletConc << endl;
  cout << "recircInletConc: "<<recircInletConc << endl;
  cout << "coreInletIndex: " << coreInletIndex << endl;
  cout << "recircInletIndex: " << recircInletIndex << endl;

  // Update variables with array splicing
  if (mats->posVelocity)
  {
    inletConc = recircConc(Eigen::seq(coreInletIndex-1,coreInletIndex),\
      Eigen::all);
    recircInletConc = dnpConc(Eigen::seq(recircInletIndex-1,recircInletIndex),\
      Eigen::all);
  } else
  {
    inletConc = recircConc(Eigen::seq(coreInletIndex,coreInletIndex+1),\
      Eigen::all);
    recircInletConc = dnpConc(Eigen::seq(recircInletIndex,recircInletIndex+1),\
      Eigen::all);
  }
  
  cout << "inletConc:"<<inletConc << endl;
  cout << "recircInletConc: "<<recircInletConc << endl;
  
  outletConc = recircConc.row(coreOutletIndex);   
  inletVelocity = mats->flowVelocity.row(recircInletIndex);
  recircOutletConc = dnpConc.row(recircOutletIndex);   
  recircInletVelocity = mats->flowVelocity.row(recircInletIndex);
  
};
//==============================================================================


