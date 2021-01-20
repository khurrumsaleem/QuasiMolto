// File: SingleGroupDNP.cpp     
// Purpose: Contains information and builds linear system for a precursor group
// Date: April 9, 2020

#include "SingleGroupDNP.h"
#include "MultiGroupDNP.h"
#include "MultiPhysicsCoupledQD.h"
#include "GreyGroupQD.h"

using namespace std;

//==============================================================================
/// SingleGroupDNP class object constructor
///
/// @param [in] myMaterials materials object for the simulation
/// @param [in] myMesh mesh object for the simulation
/// @param [in] myMGDNP multigroup dnp object this object belongs to
/// @param [in] myQD multiphysics coupled QD object for the simulation
/// @param [in] myBeta DNP fractions for this group in each neutron energy
///               group
/// @param [in] myLambda decay constant for this group 
/// @param [in] myCoreIndexOffset offset used in forming linear system in core
/// @param [in] myRecircIndexOffset offset used in forming linear system in
///               recirculation loop 
/// @param [in] myDNPid group ID for this group
SingleGroupDNP::SingleGroupDNP(Materials * myMats,\
    Mesh * myMesh,\
    MultiGroupDNP * myMGDNP,\
    Eigen::VectorXd myBeta,\
    double myLambda,\
    int myCoreIndexOffset,\
    int myRecircIndexOffset,\
    int myDNPid)
{
  // Variables to hold optional parameters
  double inpConc0;

  // Assign inputs to their member variables
  mats = myMats;
  mesh = myMesh;
  mgdnp = myMGDNP;
  beta = myBeta;
  lambda = myLambda;
  coreIndexOffset = myCoreIndexOffset;
  recircIndexOffset = myRecircIndexOffset;
  dnpID = myDNPid;

  // Initialize size of matrices and vectors
  dnpConc.setZero(mesh->nZ,mesh->nR);

  // Check for optional input  
  if ((*(mgdnp->input))["parameters"]["initial dnp concentration"])
  {

    // Read in input
    inpConc0=(*(mgdnp->input))["parameters"]["initial dnp concentration"]\
             .as<double>();

    //dnpConc = getInitialConc(inpConc0);
    dnpConc.setConstant(inpConc0);

  }

  flux.setZero(mesh->nZ+1,mesh->nR); 
  dirac.setZero(mesh->nZ+1,mesh->nR);
  recircConc.setZero(mesh->nZrecirc,mesh->nR); 

  // Check for optional input  
  if ((*(mgdnp->input))["parameters"]["initial recirc dnp concentration"])
  {

    // Read in input
    inpConc0=(*(mgdnp->input))["parameters"]\
             ["initial recirc dnp concentration"].as<double>();

    recircConc.setConstant(inpConc0);

  }

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
/// Build linear system for this precursor group. Utilized for building the core
///   and recirculation linear system. 
///
/// @param [in] myA pointer to linear system to build in
/// @param [in] myb pointer to RHS of linear system
/// @param [in] myDNPConc DNP concentration at last time step
/// @param [in] myDNPFlux DNP fluxes for modeling axial advection 
/// @param [in] myDNPFlux DNP fluxes for modeling axial advection 
/// @param [in] dzs axial heights on advecting mesh
/// @param [in] myIndexOffset row to start building linear system on 
/// @param [in] fluxSource indicator for whether a flux source is present 
void SingleGroupDNP::buildLinearSystem(Eigen::SparseMatrix<double,Eigen::RowMajor> * myA,\
    Eigen::VectorXd * myb,\
    Eigen::MatrixXd myDNPConc,\
    Eigen::MatrixXd myDNPFlux,\
    arma::rowvec dzs,\
    int myIndexOffset,\
    bool fluxSource)
{

  //int myIndex,iEq = myIndexOffset;
  int myIndex,iEq = myIndexOffset;
  int iEqTemp=0,nDNPUnknowns = myDNPConc.rows()*myDNPConc.cols();
  double coeff;
  Atemp.resize(nDNPUnknowns,myA->cols());
  Atemp.reserve(2*nDNPUnknowns);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> testMat;
  testMat.resize(nDNPUnknowns,myA->cols());
  testMat.setZero();

  #pragma omp parallel for private(myIndex,iEq,iEqTemp)
  for (int iZ = 0; iZ < myDNPConc.rows(); iZ++)
  {
    for (int iR = 0; iR < myDNPConc.cols(); iR++)
    {
      myIndex = getIndex(iZ,iR,myIndexOffset);     
      iEq = getIndex(iZ,iR,myIndexOffset);     
      iEqTemp = getIndex(iZ,iR,0);     

      //Atemp.insert(iEqTemp,myIndex) = 1 + mesh->dt*lambda; 
      testMat(iEqTemp,myIndex) = 1 + mesh->dt*lambda; 

      // Time term
      (*myb)(iEq) = myDNPConc(iZ,iR);

      // Flux source term 
      if (fluxSource)
      {
        coeff = -mesh->dt*mats->oneGroupXS->dnpFluxCoeff(iZ,iR,dnpID); 
        mgdnp->mpqd->fluxSource(iZ,iR,iEqTemp,coeff,&testMat);
        //mgdnp->mpqd->fluxSource(iZ,iR,iEqTemp,coeff,&Atemp);
      }

      // Advection term
      (*myb)(iEq) += (mesh->dt/dzs(iZ))*(myDNPFlux(iZ,iR)-myDNPFlux(iZ+1,iR));


    }
  }
  
  //myA->middleRows(myIndexOffset,nDNPUnknowns) = Atemp; 
  myA->middleRows(myIndexOffset,nDNPUnknowns) = testMat.sparseView(); 
};
//==============================================================================

//==============================================================================
/// Build linear system for this precursor group. Utilized for building the core
///   and recirculation linear system. 
///
/// @param [in] myA pointer to linear system to build in
/// @param [in] myb pointer to RHS of linear system
/// @param [in] myDNPConc DNP concentration at last time step
/// @param [in] myDNPFlux DNP fluxes for modeling axial advection 
/// @param [in] myDNPFlux DNP fluxes for modeling axial advection 
/// @param [in] dzs axial heights on advecting mesh
/// @param [in] myIndexOffset row to start building linear system on 
/// @param [in] fluxSource indicator for whether a flux source is present 
void SingleGroupDNP::buildSteadyStateLinearSystem(\
    Eigen::SparseMatrix<double,Eigen::RowMajor> * myA,\
    Eigen::VectorXd * myb,\
    Eigen::MatrixXd myDNPConc,\
    Eigen::MatrixXd myDNPFlux,\
    Eigen::MatrixXd myInletDNP,\
    arma::rowvec dzs,\
    int myIndexOffset,\
    bool fluxSource)
{

  //int myIndex,iEq = myIndexOffset;
  int upwindIndex,myIndex,iEq = myIndexOffset;
  int iEqTemp=0,nDNPUnknowns = myDNPConc.rows()*myDNPConc.cols();
  double fissionCoeff,keff,neutronFlux;
  Atemp.resize(nDNPUnknowns,myA->cols());
  Atemp.reserve(2*nDNPUnknowns);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> testMat;
  testMat.resize(nDNPUnknowns,myA->cols());
  testMat.setZero();


  if (mats->posVelocity) 
  {
    #pragma omp parallel for private(myIndex,iEq,iEqTemp)
    for (int iZ = 0; iZ < myDNPConc.rows(); iZ++)
    {
      for (int iR = 0; iR < myDNPConc.cols(); iR++)
      {
        myIndex = getIndex(iZ,iR,myIndexOffset);     
        upwindIndex = getIndex(iZ-1,iR,myIndexOffset);     
        iEq = getIndex(iZ,iR,myIndexOffset);     
        iEqTemp = getIndex(iZ,iR,0);     

        // DNP decay term
        testMat(iEqTemp,myIndex) = lambda; 

        // Flux source term 
        if (fluxSource)
        {
          fissionCoeff = mats->oneGroupXS->dnpFluxCoeff(iZ,iR,dnpID); 
          keff = mats->oneGroupXS->keff; 
          neutronFlux = mgdnp->mpqd->ggqd->sFlux(iZ,iR);
          (*myb)(iEq) = fissionCoeff*neutronFlux/keff;
        }

        // Advection terms

        // Upwind cell
        if (iZ == 0) // boundary case
          (*myb)(iEq) += myDNPFlux(iZ,iR)*myInletDNP(1,iR)/dzs(iZ);
        else
          testMat(iEqTemp,upwindIndex) = -myDNPFlux(iZ,iR)/dzs(iZ);

        // Primary cell
        testMat(iEqTemp,myIndex) += myDNPFlux(iZ+1,iR)/dzs(iZ);
      }
    }
  }
  else
  {
#pragma omp parallel for private(myIndex,iEq,iEqTemp)
    for (int iZ = 0; iZ < myDNPConc.rows(); iZ++)
    {
      for (int iR = 0; iR < myDNPConc.cols(); iR++)
      {
        myIndex = getIndex(iZ,iR,myIndexOffset);     
        upwindIndex = getIndex(iZ+1,iR,myIndexOffset);     
        iEq = getIndex(iZ,iR,myIndexOffset);     
        iEqTemp = getIndex(iZ,iR,0);     

        testMat(iEqTemp,myIndex) = lambda; 

        // Flux source term 
        if (fluxSource)
        {
          fissionCoeff = mats->oneGroupXS->dnpFluxCoeff(iZ,iR,dnpID); 
          keff = mats->oneGroupXS->keff; 
          neutronFlux = mgdnp->mpqd->ggqd->sFlux(iZ,iR);
          (*myb)(iEq) = fissionCoeff*neutronFlux/keff;
        }

        // Advection terms
    
        // Upwind cell
        if (iZ == myDNPConc.rows()) // boundary case
          (*myb)(iEq) -= myDNPFlux(iZ+1,iR)*myInletDNP(0,iR)/dzs(iZ);
        else
          testMat(iEqTemp,upwindIndex) = myDNPFlux(iZ+1,iR)/dzs(iZ);
    
        // Primary cell
        testMat(iEqTemp,myIndex) -= myDNPFlux(iZ,iR)/dzs(iZ);
      }
    }
  }

  // Splice testMat into linear system
  myA->middleRows(myIndexOffset,nDNPUnknowns) = testMat.sparseView(); 
};
//==============================================================================


//==============================================================================
/// Calculate diracs to model advection of precursors
///
/// @param [in] dnpConc DNP concentration from last time step
/// @param [in] inletConc DNP concentration at inlet
/// @param [in] outletConc DNP concentration at outlet
/// @param [out] myDirac diracs used to calculat DNP flux
Eigen::MatrixXd SingleGroupDNP::calcDiracs(Eigen::MatrixXd dnpConc,\
    Eigen::MatrixXd inletConc,\
    Eigen::VectorXd outletConc)
{

  // Declare temporary variables
  Eigen::MatrixXd myDirac(dnpConc.rows()+1,dnpConc.cols());
  double CupwindInterface,Cinterface,theta,phi;
  int lastDiracIndex = myDirac.rows()-1;

  // If flux is positive
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
  } 
  // If flux is negative
  else 
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

  return myDirac;

};
//==============================================================================

//==============================================================================
/// Calculate fluxes to model advection of precursors
///
/// @param [in] myDNPConc DNP concentration at last time step
/// @param [in] myFlowVelocity flow velocity acting on precursors  
/// @param [in] myDirac diracs to be used in flux calculation
/// @param [in] myInletConc DNP concentration at inlet
/// @param [in] myInletVelocity velocity at inlet
/// @param [in] dzs axial heights on advecting mesh
/// @param [out] myFlux fluxes used to model precursor advection
Eigen::MatrixXd SingleGroupDNP::calcFluxes(Eigen::MatrixXd myDNPConc,\
    Eigen::MatrixXd myFlowVelocity,\
    Eigen::MatrixXd myDirac,\
    Eigen::MatrixXd myInletConc,\
    Eigen::VectorXd myInletVelocity,\
    arma::rowvec dzs)
{

  // Declare temporary variables
  Eigen::MatrixXd myFlux;
  double vel; // shorthand for velocity
  int lastFluxIndex = myDirac.rows()-1;

  // Initialize DNP flux matrix
  myFlux.setZero(myDirac.rows(),myDirac.cols());

  // If velocity is positive
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

  } 
  // If velocity is negative
  else
  {

    for (int iR = 0; iR < myFlux.cols(); iR++)
    {

      // Handle all other cases
      for (int iZ = 0; iZ < myFlux.rows()-1; iZ++)
      {
        vel = myFlowVelocity(iZ,iR);
        myFlux(iZ,iR) = vel*myDNPConc(iZ,iR)\
                        + 0.5*abs(vel)*(1-abs(vel*mesh->dt/dzs(iZ)))\
                        *myDirac(iZ,iR);
      }

      // Handle iZ = nZ case
      vel = myInletVelocity(iR);
      myFlux(lastFluxIndex,iR) = vel*myInletConc(0,iR)\
                                 + 0.5*abs(vel)\
                                 *(1-abs(vel*mesh->dt/dzs(lastFluxIndex-1)))\
                                 *myDirac(lastFluxIndex,iR);

    }
  }

  return myFlux;

};
//==============================================================================

//==============================================================================
/// Calculate implicit fluxes to model advection of precursors
///
/// @param [in] myDNPConc DNP concentration at last time step
/// @param [in] myFlowVelocity flow velocity acting on precursors  
/// @param [in] myDirac diracs to be used in flux calculation
/// @param [in] myInletConc DNP concentration at inlet
/// @param [in] myInletVelocity velocity at inlet
/// @param [in] dzs axial heights on advecting mesh
/// @param [out] myFlux fluxes used to model precursor advection
Eigen::MatrixXd SingleGroupDNP::calcImplicitFluxes(Eigen::MatrixXd myDNPConc,\
    Eigen::MatrixXd myFlowVelocity,\
    Eigen::MatrixXd myInletConc,\
    Eigen::VectorXd myInletVelocity,
    arma::rowvec dzs)
{

  // Declare temporary variables
  Eigen::MatrixXd myFlux(myDNPConc.rows()+1,myDNPConc.cols());
  double vel; // shorthand for velocity
  int lastFluxIndex = myFlux.rows()-1;

  // If velocity is positive
  if (mats->posVelocity) {

    for (int iR = 0; iR < myFlux.cols(); iR++)
    {
      // Handle iZ = 0 case
      vel = myInletVelocity(iR);
      myFlux(0,iR) = vel;

      // Handle all other cases
      for (int iZ = 1; iZ < myFlux.rows(); iZ++)
      {
        vel = myFlowVelocity(iZ-1,iR);
        myFlux(iZ,iR) = vel;
      }

    }

  } 
  // If velocity is negative
  else
  {

    for (int iR = 0; iR < myFlux.cols(); iR++)
    {

      // Handle all other cases
      for (int iZ = 0; iZ < myFlux.rows()-1; iZ++)
      {
        vel = myFlowVelocity(iZ,iR);
        myFlux(iZ,iR) = vel;
      }

      // Handle iZ = nZ case
      vel = myInletVelocity(iR);
      myFlux(lastFluxIndex,iR) = vel;

    }
  }

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
int SingleGroupDNP::getCoreConc()
{
  PetscErrorCode ierr;
  PetscScalar value[1]; 
  PetscInt index[1]; 
  VecScatter     ctx;
  Vec temp_x_p_seq;

  if (mesh->petsc)
  {

    // Initialize temporary vector
    initPETScVec(&temp_x_p_seq,mgdnp->mpqd->nUnknowns);

    // Gather values of x_p on all procs
    VecScatterCreateToAll(mgdnp->mpqd->x_p,&ctx,&temp_x_p_seq);
    VecScatterBegin(ctx,mgdnp->mpqd->x_p,temp_x_p_seq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,mgdnp->mpqd->x_p,temp_x_p_seq,INSERT_VALUES,SCATTER_FORWARD);

    // loop over spatial mesh
    for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
    {
      for (int iZ = 0; iZ < mesh->dzsCorner.size(); iZ++)
      {

        // Read fluxes into flux vector
        index[0] = getIndex(iZ,iR,coreIndexOffset);
        ierr = VecGetValues(temp_x_p_seq,1,index,value);
        dnpConc(iZ,iR) = value[0];

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

        dnpConc(iZ,iR) = mgdnp->mpqd->x(getIndex(iZ,iR,coreIndexOffset));   

      }
    }
  }

};
//==============================================================================

//==============================================================================
/// Map solution from 2D matrix to 1D vector
///
void SingleGroupDNP::setCoreConc()
{

  for (int iR = 0; iR < mesh-> drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh-> dzsCorner.size(); iZ++)
    { 

      mgdnp->mpqd->xPast(getIndex(iZ,iR,coreIndexOffset)) = dnpConc(iZ,iR);   

    }
  }

};
//==============================================================================

//==============================================================================
/// Map solution in 1D vector to 2D solution
///
int SingleGroupDNP::getRecircConc()
{

  PetscErrorCode ierr;
  PetscScalar value[1]; 
  PetscInt index[1]; 
  VecScatter     ctx;
  Vec temp_x_p_seq;

  if (mesh->petsc)
  {

    // Initialize temporary vector
    initPETScVec(&temp_x_p_seq,mgdnp->mpqd->nUnknowns);

    // Gather values of x_p on all procs
    VecScatterCreateToAll(mgdnp->recircx_p,&ctx,&temp_x_p_seq);
    VecScatterBegin(ctx,mgdnp->recircx_p,temp_x_p_seq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,mgdnp->recircx_p,temp_x_p_seq,INSERT_VALUES,SCATTER_FORWARD);

    // loop over spatial mesh
    for (int iR = 0; iR < mesh->drsCorner.size(); iR++)
    {
      for (int iZ = 0; iZ < mesh->nZrecirc; iZ++)
      {

        // Read fluxes into flux vector
        index[0] = getIndex(iZ,iR,recircIndexOffset);
        ierr = VecGetValues(temp_x_p_seq,1,index,value);
        recircConc(iZ,iR) = value[0];

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
      for (int iZ = 0; iZ < mesh-> nZrecirc; iZ++)
      { 

        recircConc(iZ,iR) = mgdnp->recircx(getIndex(iZ,iR,recircIndexOffset));   

      }
    }
  }

};
//==============================================================================

//==============================================================================
/// Map solution to 1D vector from 2D solution
///
void SingleGroupDNP::setRecircConc()
{

  for (int iR = 0; iR < mesh-> drsCorner.size(); iR++)
  {
    for (int iZ = 0; iZ < mesh-> nZrecirc; iZ++)
    { 

      mgdnp->recircx(getIndex(iZ,iR,recircIndexOffset)) = recircConc(iZ,iR);   

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
/// Build linear system governing transient core DNP concentrations
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
      mesh->dzsCorner,\
      coreIndexOffset);
};
//==============================================================================

//==============================================================================
/// Build linear system governing steady state core DNP concentrations
///
void SingleGroupDNP::buildSteadyStateCoreLinearSystem()
{

  Eigen::MatrixXd coreFlux;

  updateBoundaryConditions();

  coreFlux = calcImplicitFluxes(dnpConc,\
      mats->flowVelocity,\
      inletConc,\
      inletVelocity,\
      mesh->dzsCorner);

  buildSteadyStateLinearSystem(&(mgdnp->mpqd->A),\
      &(mgdnp->mpqd->b),\
      dnpConc,\
      coreFlux,\
      inletConc,\
      mesh->dzsCorner,\
      coreIndexOffset);
};
//==============================================================================

//==============================================================================
/// Build linear system governing transient recirulation DNP concentrations
///
void SingleGroupDNP::buildRecircLinearSystem()
{

  Eigen::MatrixXd recircDirac,recircFlux,dumbySigF;

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
      mesh->dzsCornerRecirc,\
      recircIndexOffset,\
      false);

};
//==============================================================================

//==============================================================================
/// Build linear system governing steady state recirulation DNP concentrations
///
void SingleGroupDNP::buildSteadyStateRecircLinearSystem()
{

  Eigen::MatrixXd recircFlux,dumbySigF;

  updateBoundaryConditions();

  recircFlux = calcImplicitFluxes(recircConc,\
      mats->recircFlowVelocity,\
      recircInletConc,\
      recircInletVelocity,\
      mesh->dzsCornerRecirc);

  buildSteadyStateLinearSystem(&(mgdnp->recircA),\
      &(mgdnp->recircb),\
      recircConc,\
      recircFlux,\
      recircInletConc,\
      mesh->dzsCornerRecirc,\
      recircIndexOffset,\
      false);

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
Eigen::MatrixXd SingleGroupDNP::getInitialConc(double initConc)
{

  Eigen::MatrixXd conc; 

  conc.setZero(mesh->nZ,mesh->nR);    

  // Check if a material is stationary. If so, assume it has no initial 
  // precursor concentration.  
  for (int iR = 0; iR < mesh->nR; iR++)
  {
    for (int iZ = 0; iZ < mesh->nZ; iZ++)
    {
      //if(not mats->matBank[mats->matMap(iZ,iR)]->stationary) 
      // conc(iZ,iR) = initConc; 
    }
  }

  return conc;
};
//==============================================================================

//==============================================================================
/// Update boundary conditions 
///
void SingleGroupDNP::updateBoundaryConditions()
{

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

  outletConc = recircConc.row(coreOutletIndex);   
  inletVelocity = mats->flowVelocity.row(recircInletIndex);
  recircOutletConc = dnpConc.row(recircOutletIndex);   
  recircInletVelocity = mats->flowVelocity.row(recircInletIndex);

};
//==============================================================================

/* PETSc */

/* STEADY STATE */

//==============================================================================
/// Build linear system governing steady state core DNP concentrations
///
void SingleGroupDNP::buildSteadyStateCoreLinearSystem_p()
{

  Eigen::MatrixXd coreFlux;

  updateBoundaryConditions();

  coreFlux = calcImplicitFluxes(dnpConc,\
      mats->flowVelocity,\
      inletConc,\
      inletVelocity,\
      mesh->dzsCorner);

  buildSteadyStateLinearSystem_p(&(mgdnp->mpqd->A_p),\
      &(mgdnp->mpqd->b_p),\
      dnpConc,\
      coreFlux,\
      inletConc,\
      mesh->dzsCorner,\
      coreIndexOffset);
};
//==============================================================================

//==============================================================================
/// Build linear system governing steady state recirulation DNP concentrations
///
void SingleGroupDNP::buildSteadyStateRecircLinearSystem_p()
{

  Eigen::MatrixXd recircFlux,dumbySigF;

  updateBoundaryConditions();

  recircFlux = calcImplicitFluxes(recircConc,\
      mats->recircFlowVelocity,\
      recircInletConc,\
      recircInletVelocity,\
      mesh->dzsCornerRecirc);

  buildSteadyStateLinearSystem_p(&(mgdnp->recircA_p),\
      &(mgdnp->recircb_p),\
      recircConc,\
      recircFlux,\
      recircInletConc,\
      mesh->dzsCornerRecirc,\
      recircIndexOffset,\
      false);

};
//==============================================================================

//==============================================================================
/// Build linear system for this precursor group. Utilized for building the core
///   and recirculation linear system. 
///
/// @param [in] myA pointer to linear system to build in
/// @param [in] myb pointer to RHS of linear system
/// @param [in] myDNPConc DNP concentration at last time step
/// @param [in] myDNPFlux DNP fluxes for modeling axial advection 
/// @param [in] myDNPFlux DNP fluxes for modeling axial advection 
/// @param [in] dzs axial heights on advecting mesh
/// @param [in] myIndexOffset row to start building linear system on 
/// @param [in] fluxSource indicator for whether a flux source is present 
int SingleGroupDNP::buildSteadyStateLinearSystem_p(\
    Mat * A_p,\
    Vec * b_p,\
    Eigen::MatrixXd myDNPConc,\
    Eigen::MatrixXd myDNPFlux,\
    Eigen::MatrixXd myInletDNP,\
    arma::rowvec dzs,\
    int myIndexOffset,\
    bool fluxSource)
{

  //int myIndex,iEq = myIndexOffset;
  int upwindIndex,myIndex,iEq = myIndexOffset;
  int iEqTemp=0,nDNPUnknowns = myDNPConc.rows()*myDNPConc.cols();
  double fissionCoeff,keff,neutronFlux;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> testMat;
  PetscErrorCode ierr;
  PetscScalar value;

  //testMat.resize(nDNPUnknowns,myA->cols());
  //testMat.setZero();
  //Atemp.resize(nDNPUnknowns,myA->cols());
  //Atemp.reserve(2*nDNPUnknowns);

  if (mats->posVelocity) 
  {
    ////    #pragma omp parallel for private(myIndex,iEq,iEqTemp)
    for (int iZ = 0; iZ < myDNPConc.rows(); iZ++)
    {
      for (int iR = 0; iR < myDNPConc.cols(); iR++)
      {
        myIndex = getIndex(iZ,iR,myIndexOffset);     
        upwindIndex = getIndex(iZ-1,iR,myIndexOffset);     
        iEq = getIndex(iZ,iR,myIndexOffset);     
        iEqTemp = getIndex(iZ,iR,0);     

        // DNP decay term
        ierr = MatSetValue(*A_p,iEq,myIndex,lambda,ADD_VALUES);CHKERRQ(ierr); 
        //testMat(iEqTemp,myIndex) = lambda; 

        // Flux source term 
        if (fluxSource)
        {
          fissionCoeff = mats->oneGroupXS->dnpFluxCoeff(iZ,iR,dnpID); 
          keff = mats->oneGroupXS->keff; 
          neutronFlux = mgdnp->mpqd->ggqd->sFlux(iZ,iR);
          value = fissionCoeff*neutronFlux/keff;
          ierr = VecSetValue(*b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
          //(*myb)(iEq) = fissionCoeff*neutronFlux/keff;
        }

        // Advection terms

        // Upwind cell
        if (iZ == 0) // boundary case
        {
          value = myDNPFlux(iZ,iR)*myInletDNP(1,iR)/dzs(iZ);
          ierr = VecSetValue(*b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
          //(*myb)(iEq) += myDNPFlux(iZ,iR)*myInletDNP(1,iR)/dzs(iZ);
        }
        else
        {
          value = -myDNPFlux(iZ,iR)/dzs(iZ);
          ierr = MatSetValue(*A_p,iEq,upwindIndex,value,ADD_VALUES);CHKERRQ(ierr); 
          //testMat(iEqTemp,upwindIndex) = -myDNPFlux(iZ,iR)/dzs(iZ);
        }

        // Primary cell
        value = myDNPFlux(iZ+1,iR)/dzs(iZ);
        ierr = MatSetValue(*A_p,iEq,myIndex,value,ADD_VALUES);CHKERRQ(ierr); 
        //testMat(iEqTemp,myIndex) += myDNPFlux(iZ+1,iR)/dzs(iZ);
      }
    }
  }
  else
  {
    //#pragma omp parallel for private(myIndex,iEq,iEqTemp)
    for (int iZ = 0; iZ < myDNPConc.rows(); iZ++)
    {
      for (int iR = 0; iR < myDNPConc.cols(); iR++)
      {
        myIndex = getIndex(iZ,iR,myIndexOffset);     
        upwindIndex = getIndex(iZ+1,iR,myIndexOffset);     
        iEq = getIndex(iZ,iR,myIndexOffset);     
        iEqTemp = getIndex(iZ,iR,0);     

        //testMat(iEqTemp,myIndex) = lambda; 
        ierr = MatSetValue(*A_p,iEq,myIndex,lambda,ADD_VALUES);CHKERRQ(ierr); 

        // Flux source term 
        if (fluxSource)
        {
          fissionCoeff = mats->oneGroupXS->dnpFluxCoeff(iZ,iR,dnpID); 
          keff = mats->oneGroupXS->keff; 
          neutronFlux = mgdnp->mpqd->ggqd->sFlux(iZ,iR);
          value = fissionCoeff*neutronFlux/keff;
          ierr = VecSetValue(*b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
          //(*myb)(iEq) = fissionCoeff*neutronFlux/keff;
        }

        // Advection terms

        // Upwind cell
        if (iZ == myDNPConc.rows()) // boundary case
        {
          value = -myDNPFlux(iZ+1,iR)*myInletDNP(0,iR)/dzs(iZ);
          ierr = VecSetValue(*b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
          //(*myb)(iEq) -= myDNPFlux(iZ+1,iR)*myInletDNP(0,iR)/dzs(iZ);
        }
        else
        {
          value = myDNPFlux(iZ+1,iR)/dzs(iZ);
          ierr = MatSetValue(*A_p,iEq,upwindIndex,value,ADD_VALUES);CHKERRQ(ierr); 
          //testMat(iEqTemp,upwindIndex) = myDNPFlux(iZ+1,iR)/dzs(iZ);
        }

        // Primary cell
        value = -myDNPFlux(iZ,iR)/dzs(iZ);
        ierr = MatSetValue(*A_p,iEq,myIndex,value,ADD_VALUES);CHKERRQ(ierr); 
        //testMat(iEqTemp,myIndex) -= myDNPFlux(iZ,iR)/dzs(iZ);
      }
    }
  }

  // Splice testMat into linear system
  //myA->middleRows(myIndexOffset,nDNPUnknowns) = testMat.sparseView(); 
};
//==============================================================================

/* TRANSIENT */

//==============================================================================
/// Build linear system governing transient core DNP concentrations
///
void SingleGroupDNP::buildCoreLinearSystem_p()
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

  buildLinearSystem_p(&(mgdnp->mpqd->A_p),\
      &(mgdnp->mpqd->b_p),\
      dnpConc,\
      coreFlux,\
      mesh->dzsCorner,\
      coreIndexOffset);
};
//==============================================================================

//==============================================================================
/// Build linear system governing transient recirulation DNP concentrations
///
void SingleGroupDNP::buildRecircLinearSystem_p()
{

  Eigen::MatrixXd recircDirac,recircFlux,dumbySigF;

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

  buildLinearSystem_p(&(mgdnp->recircA_p),\
      &(mgdnp->recircb_p),\
      recircConc,\
      recircFlux,\
      mesh->dzsCornerRecirc,\
      recircIndexOffset,\
      false);

};
//==============================================================================

//==============================================================================
/// Build linear system for this precursor group. Utilized for building the core
///   and recirculation linear system. 
///
/// @param [in] myA pointer to linear system to build in
/// @param [in] myb pointer to RHS of linear system
/// @param [in] myDNPConc DNP concentration at last time step
/// @param [in] myDNPFlux DNP fluxes for modeling axial advection 
/// @param [in] myDNPFlux DNP fluxes for modeling axial advection 
/// @param [in] dzs axial heights on advecting mesh
/// @param [in] myIndexOffset row to start building linear system on 
/// @param [in] fluxSource indicator for whether a flux source is present 
int SingleGroupDNP::buildLinearSystem_p(
    Mat * A_p,\
    Vec * b_p,\
    Eigen::MatrixXd myDNPConc,\
    Eigen::MatrixXd myDNPFlux,\
    arma::rowvec dzs,\
    int myIndexOffset,\
    bool fluxSource)
{

  //int myIndex,iEq = myIndexOffset;
  int myIndex,iEq = myIndexOffset;
  int iEqTemp=0,nDNPUnknowns = myDNPConc.rows()*myDNPConc.cols();
  double coeff;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> testMat;
  PetscErrorCode ierr;
  PetscScalar value;
  
  //#pragma omp parallel for private(myIndex,iEq,iEqTemp)
  for (int iZ = 0; iZ < myDNPConc.rows(); iZ++)
  {
    for (int iR = 0; iR < myDNPConc.cols(); iR++)
    {
      myIndex = getIndex(iZ,iR,myIndexOffset);     
      iEq = getIndex(iZ,iR,myIndexOffset);     
      iEqTemp = getIndex(iZ,iR,0);     

      value = 1 + mesh->dt*lambda;
      ierr = MatSetValue(*A_p,iEq,myIndex,lambda,ADD_VALUES);CHKERRQ(ierr); 
      //testMat(iEqTemp,myIndex) = 1 + mesh->dt*lambda; 

      // Time term
      value = myDNPConc(iZ,iR);
      ierr = VecSetValue(*b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
      //(*myb)(iEq) = myDNPConc(iZ,iR);

      // Flux source term 
      if (fluxSource)
      {
        coeff = -mesh->dt*mats->oneGroupXS->dnpFluxCoeff(iZ,iR,dnpID); 
        mgdnp->mpqd->fluxSource(iZ,iR,iEqTemp,coeff,&testMat);
      }

      // Advection term
      value = (mesh->dt/dzs(iZ))*(myDNPFlux(iZ,iR)-myDNPFlux(iZ+1,iR));
      ierr = VecSetValue(*b_p,iEq,value,ADD_VALUES);CHKERRQ(ierr); 
      //(*myb)(iEq) += (mesh->dt/dzs(iZ))*(myDNPFlux(iZ,iR)-myDNPFlux(iZ+1,iR));


    }
  }
  
  //myA->middleRows(myIndexOffset,nDNPUnknowns) = Atemp; 
  //myA->middleRows(myIndexOffset,nDNPUnknowns) = testMat.sparseView(); 
};
//==============================================================================

