// File: TransportToQDCoupling.cpp
// Purpose: handle communication between transport and quasidiffusion objects
// Date: February 12, 2020

#include "TransportToQDCoupling.h"
#include "SimpleCornerBalance.h"

using namespace std;

//==============================================================================
/// TransportToQDCoupling class object constructor
///
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myInput YAML input object for the simulation
/// @param [in] myMGT Multigroup transport object for the simulation
/// @param [in] myMGQD Multigroup quasidiffusion object for the simulation
TransportToQDCoupling::TransportToQDCoupling(Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput,\
    MultiGroupTransport * myMGT,\
    MultiGroupQD * myMGQD)
{
  // Assign pointers for materials, mesh, input objects, multigroup transport
  // and quasidiffusion objects
  materials = myMaterials;
  mesh = myMesh;
  input = myInput;
  MGT = myMGT;
  MGQD = myMGQD;

  // Check for optional parameters
  checkOptionalParams(); 

};
//==============================================================================

//==============================================================================
/// Calculate Eddington factors using angular fluxes from transport objects
/// @param [out] allConverged boolean indicating if the residual on the 
/// Eddington factors passes the convergence criteria
bool TransportToQDCoupling::calcEddingtonFactors()
{
  int rows = MGT->SGTs[0]->sFlux.rows();
  int cols = MGT->SGTs[0]->sFlux.cols();
  int angIdx,xiIdx=0,muIdx=1,etaIdx=2,weightIdx = 3;  
  double angFlux,mu,xi,weight,EzzCoef,ErrCoef,ErzCoef;
  double numeratorEzz,numeratorErr,numeratorErz,denominator;
  double residualZz,residualRr,residualRz;
  bool interfaceConverged,cellAvgConverged=true;

  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    // store past eddington factors
    MGQD->SGQDs[iGroup]->EzzPrev = MGQD->SGQDs[iGroup]->Ezz;
    MGQD->SGQDs[iGroup]->ErrPrev = MGQD->SGQDs[iGroup]->Err;
    MGQD->SGQDs[iGroup]->ErzPrev = MGQD->SGQDs[iGroup]->Erz;

    for (int iR = 0; iR < cols; iR++)
    {
      for (int iZ = 0; iZ < rows; iZ++)
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

    // measure the residual for each Eddington factors 
    residualZz = ((MGQD->SGQDs[iGroup]->Ezz - MGQD->SGQDs[iGroup]->EzzPrev)\
        .cwiseQuotient(MGQD->SGQDs[iGroup]->Ezz)).norm();
    residualRr = ((MGQD->SGQDs[iGroup]->Err - MGQD->SGQDs[iGroup]->ErrPrev)\
        .cwiseQuotient(MGQD->SGQDs[iGroup]->Err)).norm();
    residualRz = ((MGQD->SGQDs[iGroup]->Erz - MGQD->SGQDs[iGroup]->ErzPrev)\
        .cwiseQuotient(MGQD->SGQDs[iGroup]->Erz)).norm();

    residualZz = calcResidual(MGQD->SGQDs[iGroup]->EzzPrev,\
        MGQD->SGQDs[iGroup]->Ezz);
    residualRr = calcResidual(MGQD->SGQDs[iGroup]->ErrPrev,\
        MGQD->SGQDs[iGroup]->Err);
    //residualRz = calcResidual(MGQD->SGQDs[iGroup]->ErzPrev,\
        MGQD->SGQDs[iGroup]->Erz);

   // cout << "residualZz: " << residualZz << endl;
   // cout << endl;
   // cout << "residualRr: " << residualRr << endl;
   // cout << endl;
   // cout << "residualRz: " << residualRz << endl;
   // cout << endl;

    //if (residualZz < epsEddington and residualRr < epsEddington and 
    //    residualRz < epsEddington)
    if (residualZz < epsEddington and residualRr < epsEddington)
      cellAvgConverged = cellAvgConverged and true;
    else
      cellAvgConverged = cellAvgConverged and false;

  } //iGroup


  // Calculate interface Eddingtons
  interfaceConverged = calcInterfaceEddingtonFactors(); 

  // Calculate G used in integrating factor
  calcGFactors();
  
  // Calculate g0 and g1 coefficients used in integrating factor 
  calcIntFactorCoeffs();


  return (cellAvgConverged and interfaceConverged);
}
//==============================================================================

//==============================================================================
/// Calculate interfaceEddington factors using angular fluxes from transport 
///     objects
/// @param [out] allConverged boolean indicating if the residual on the 
/// Eddington factors passes the convergence criteria
bool TransportToQDCoupling::calcInterfaceEddingtonFactors()
{
  int rows = MGT->SGTs[0]->sFlux.rows();
  int cols = MGT->SGTs[0]->sFlux.cols();
  int angIdx,xiIdx=0,muIdx=1,etaIdx=2,weightIdx = 3;  
  double angFlux,mu,xi,weight,EzzCoef,ErrCoef,ErzCoef;
  double numeratorEzz,numeratorErr,numeratorErz,denominator;
  double residualZz,residualRr,residualRz;
  double volLeft,volRight,volUp,volDown;
  bool allConverged=true;

  Eigen::MatrixXd ErzAxialPrev,ErrAxialPrev,EzzAxialPrev; 
  Eigen::MatrixXd ErzRadialPrev,ErrRadialPrev,EzzRadialPrev;

  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    
    // store past eddington factors
    EzzRadialPrev = MGQD->SGQDs[iGroup]->EzzRadial;
    ErzRadialPrev = MGQD->SGQDs[iGroup]->ErzRadial;
    ErrRadialPrev = MGQD->SGQDs[iGroup]->ErrRadial;

    for (int iR = 0; iR < cols+1; iR++)
    {
      for (int iZ = 0; iZ < rows; iZ++)
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
            if (iR == 0) 
              angFlux = MGT->SGTs[iGroup]->aFlux(iZ,iR,angIdx);
            else if (iR == cols) 
              angFlux = MGT->SGTs[iGroup]->aFlux(iZ,iR-1,angIdx);
            else
            {
              volLeft = mesh->getGeoParams(iR-1,iZ)[0]; 
              volRight = mesh->getGeoParams(iR,iZ)[0];
              angFlux = (volLeft*MGT->SGTs[iGroup]->aFlux(iZ,iR-1,angIdx)\
                  + volRight*MGT->SGTs[iGroup]->aFlux(iZ,iR,angIdx))\
                        /(volLeft+volRight);
            }

            EzzCoef = xi*xi;
            ErrCoef = mu*mu;
            ErzCoef = mu*xi;

            numeratorEzz = numeratorEzz + EzzCoef*angFlux*weight;
            numeratorErr = numeratorErr + ErrCoef*angFlux*weight;
            numeratorErz = numeratorErz + ErzCoef*angFlux*weight;

            denominator = denominator + angFlux*weight;

          } //iMu
        } //iXi 

        MGQD->SGQDs[iGroup]->EzzRadial(iZ,iR) = numeratorEzz/denominator;
        MGQD->SGQDs[iGroup]->ErrRadial(iZ,iR) = numeratorErr/denominator;
        MGQD->SGQDs[iGroup]->ErzRadial(iZ,iR) = numeratorErz/denominator;

      } //iZ
    } //iR

    residualZz = calcResidual(EzzRadialPrev, MGQD->SGQDs[iGroup]->EzzRadial);
    residualRr = calcResidual(ErrRadialPrev, MGQD->SGQDs[iGroup]->ErrRadial);
    //residualRz = calcResidual(ErzRadialPrev, MGQD->SGQDs[iGroup]->ErzRadial);

   // cout << "residualZzRadial: " << residualZz << endl;
   // cout << endl;
   // cout << "residualRrRadial: " << residualRr << endl;
   // cout << endl;
   // cout << "residualRzRadial: " << residualRz << endl;
   // cout << endl;

    //if (residualZz < epsEddington and residualRr < epsEddington and 
    //    residualRz < epsEddington)
    if (residualZz < epsEddington and residualRr < epsEddington)
      allConverged = allConverged and true;
    else
      allConverged = allConverged and false;
  } //iGroup

  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    
    // store past eddington factors
    EzzAxialPrev = MGQD->SGQDs[iGroup]->EzzAxial;
    ErzAxialPrev = MGQD->SGQDs[iGroup]->ErzAxial;
    ErrAxialPrev = MGQD->SGQDs[iGroup]->ErrAxial;

    for (int iR = 0; iR < cols; iR++)
    {
      for (int iZ = 0; iZ < rows+1; iZ++)
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
            if (iZ == 0) 
              angFlux = MGT->SGTs[iGroup]->aFlux(iZ,iR,angIdx);
            else if (iZ == rows) 
              angFlux = MGT->SGTs[iGroup]->aFlux(iZ-1,iR,angIdx);
            else
            {
              volDown = mesh->getGeoParams(iR,iZ-1)[0]; 
              volUp = mesh->getGeoParams(iR,iZ)[0];
              angFlux = (volDown*MGT->SGTs[iGroup]->aFlux(iZ-1,iR,angIdx)\
                  + volUp*MGT->SGTs[iGroup]->aFlux(iZ,iR,angIdx))\
                        /(volUp+volDown);
            }

            EzzCoef = xi*xi;
            ErrCoef = mu*mu;
            ErzCoef = mu*xi;

            numeratorEzz = numeratorEzz + EzzCoef*angFlux*weight;
            numeratorErr = numeratorErr + ErrCoef*angFlux*weight;
            numeratorErz = numeratorErz + ErzCoef*angFlux*weight;

            denominator = denominator + angFlux*weight;

          } //iMu
        } //iXi 

        MGQD->SGQDs[iGroup]->EzzAxial(iZ,iR) = numeratorEzz/denominator;
        MGQD->SGQDs[iGroup]->ErrAxial(iZ,iR) = numeratorErr/denominator;
        MGQD->SGQDs[iGroup]->ErzAxial(iZ,iR) = numeratorErz/denominator;

      } //iZ
    } //iR

    // measure the residual for each Eddington factors 
    residualZz = calcResidual(EzzAxialPrev, MGQD->SGQDs[iGroup]->EzzAxial);
    residualRr = calcResidual(ErrAxialPrev, MGQD->SGQDs[iGroup]->ErrAxial);
    //residualRz = calcResidual(ErzAxialPrev, MGQD->SGQDs[iGroup]->ErzAxial);

   // cout << "residualZzAxial: " << residualZz << endl;
   // cout << endl;
   // cout << "residualRrAxial: " << residualRr << endl;
   // cout << endl;
   // cout << "residualRzAxial: " << residualRz << endl;
   // cout << endl;

    //if (residualZz < epsEddington and residualRr < epsEddington and 
    //    residualRz < epsEddington)
    if (residualZz < epsEddington and residualRr < epsEddington)
      allConverged = allConverged and true;
    else
      allConverged = allConverged and false;
  } //iGroup


  return allConverged;
}
//==============================================================================

//==============================================================================
/// Calculate G factors using group Eddington factors 
void TransportToQDCoupling::calcGFactors()
{
  int rows = MGT->SGTs[0]->sFlux.rows();
  int cols = MGT->SGTs[0]->sFlux.cols();
  double Ezz, Err;

  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    for (int iZ = 0; iZ < rows; iZ++)
    {
      for (int iR = 0; iR < cols; iR++)
      {

        // Calculate cell center G
        Ezz = MGQD->SGQDs[iGroup]->Ezz(iZ,iR);
        Err = MGQD->SGQDs[iGroup]->Err(iZ,iR);

        MGQD->SGQDs[iGroup]->G(iZ,iR)= 1.0 + (Err + Ezz - 1.0) / Err;

        // Calculate cell edge G
        Ezz = MGQD->SGQDs[iGroup]->EzzRadial(iZ,iR);
        Err = MGQD->SGQDs[iGroup]->ErrRadial(iZ,iR);

        MGQD->SGQDs[iGroup]->GRadial(iZ,iR)= 1.0 + (Err + Ezz - 1.0) / Err;

      } // iR

      // Calculate cell edge case for edge G
      Ezz = MGQD->SGQDs[iGroup]->EzzRadial(iZ,cols);
      Err = MGQD->SGQDs[iGroup]->ErrRadial(iZ,cols);

      MGQD->SGQDs[iGroup]->GRadial(iZ,cols)= 1.0 + (Err + Ezz - 1.0) / Err;

    } // iZ

  } //iGroup

}
//==============================================================================

//==============================================================================
/// Calculate integrating factor coefficients
void TransportToQDCoupling::calcIntFactorCoeffs()
{
  int rows = MGT->SGTs[0]->sFlux.rows();
  double g0, g1, Gcell, Gedge, p;

  double rAvg = mesh->rVWCornerCent[0], rUp = mesh->rCornerEdge[1];

  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    for (int iZ = 0; iZ < rows; iZ++)
    {

      // Calculate g0 and g1 coefficients and store them in SGQD objects
      Gcell = MGQD->SGQDs[iGroup]->G(iZ,0); 
      Gedge = MGQD->SGQDs[iGroup]->G(iZ,0);
      p = 2;
      g1 = (1.0 / (rUp-rAvg)) * ((Gedge/pow(rUp,p)) - (Gcell/pow(rAvg,p)));
      g0 = (Gcell/pow(rAvg,p)) - g1*rAvg;

      MGQD->SGQDs[iGroup]->g1(iZ) = g1;
      MGQD->SGQDs[iGroup]->g0(iZ) = g0;

    } // iZ
      
  } //iGroup

}
//==============================================================================

//==============================================================================
/// Calculate a number a parameters used for forming the boundary conditions of
/// low order problem 
void TransportToQDCoupling::calcBCs()
{
  int rows = MGT->SGTs[0]->sFlux.rows();
  int cols = MGT->SGTs[0]->sFlux.cols();
  int eIdx = cols - 1,sIdx = rows - 1;
  int angIdx,xiIdx=0,muIdx=1,etaIdx=2,weightIdx = 3;  
  double angFlux,angFluxN,angFluxS,mu,xi,weight;
  double inwardJrE,inwardJzN,inwardJzS;
  double inwardFluxE,inwardFluxN,inwardFluxS;
  double outwardJrE,outwardJzN,outwardJzS;
  double outwardFluxE,outwardFluxN,outwardFluxS;
  double localScalarFluxE, localScalarFluxN, localScalarFluxS;

  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    for (int iZ = 0; iZ < rows; iZ++)
    {
      // reset accumulators
      inwardJrE = 0.0;
      inwardFluxE = 0.0;
      outwardJrE = 0.0;
      outwardFluxE = 0.0;
      localScalarFluxE = 0.0;

      // loop over quadrature
      for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi)
      {
        xi = mesh->quadrature[iXi].quad[0][xiIdx];
        for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu)
        {
          angIdx=mesh->quadrature[iXi].ordIdx[iMu];
          mu = mesh->quadrature[iXi].quad[iMu][muIdx]; 
          weight = mesh->quadrature[iXi].quad[iMu][weightIdx];

          // Get angular fluxes. If at a location where there is a boundary
          // condition for the angular flux, use that. If not, grab the 
          // cell-average angular flux of the nearest cell. 
          if (mu < 0)
          {
            angFlux = MGT->SCBSolve->outerBC[iGroup];           
          }
          else
          {
            angFlux = MGT->SGTs[iGroup]->aFlux(iZ,eIdx,angIdx);
          } 

          localScalarFluxE += angFlux*weight;

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

      // set eastern flux bc
      MGQD->SGQDs[iGroup]->eFluxBC(iZ) = localScalarFluxE; 
      MGQD->SGQDs[iGroup]->eCurrentRBC(iZ) = outwardJrE+inwardJrE;

      // set absolute current
      MGQD->SGQDs[iGroup]->eAbsCurrentBC(iZ) = outwardJrE - inwardJrE;

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
      localScalarFluxN = 0.0;
      localScalarFluxS = 0.0;

      // loop over quadrature
      for (int iXi = 0; iXi < mesh->quadrature.size(); ++iXi)
      {
        xi = mesh->quadrature[iXi].quad[0][xiIdx];
        for (int iMu = 0; iMu < mesh->quadrature[iXi].nOrd; ++iMu)
        {
          angIdx=mesh->quadrature[iXi].ordIdx[iMu];
          mu = mesh->quadrature[iXi].quad[iMu][muIdx]; 
          weight = mesh->quadrature[iXi].quad[iMu][weightIdx];

          // Get angular fluxes. If at a location where there is a boundary
          // condition for the angular flux, use that. If not, grab the 
          // cell-average angular flux of the nearest cell. 
          if (xi > 0)
          {
            angFluxN = MGT->SCBSolve->lowerBC[iGroup];           
            angFluxS = MGT->SGTs[iGroup]->aFlux(sIdx,iR,angIdx);
          }
          else
          {
            angFluxN = MGT->SGTs[iGroup]->aFlux(0,iR,angIdx);
            angFluxS = MGT->SCBSolve->upperBC[iGroup];           
          } 

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

      // set flux BCs
      MGQD->SGQDs[iGroup]->nFluxBC(iR) = localScalarFluxN;
      MGQD->SGQDs[iGroup]->sFluxBC(iR) = localScalarFluxS;
      MGQD->SGQDs[iGroup]->nCurrentZBC(iR) = outwardJzN+inwardJzN;
      MGQD->SGQDs[iGroup]->sCurrentZBC(iR) = outwardJzS+inwardJzS;

      // set absolute currents
      MGQD->SGQDs[iGroup]->nAbsCurrentBC(iR) = outwardJzN - inwardJzN;
      MGQD->SGQDs[iGroup]->sAbsCurrentBC(iR) = outwardJzS - inwardJzS;

    } //iR 

  } //iGroup

}
//==============================================================================

//==============================================================================
/// Calculate residual between two matrices 
///
double TransportToQDCoupling::calcResidual(Eigen::MatrixXd matrix1,\
    Eigen::MatrixXd matrix2)
{

  Eigen::MatrixXd ones,residualMatrix;
  double residual;
  double min = 1E-12;

  ones.setOnes(matrix1.rows(),matrix1.cols());

  for (int iRow = 0; iRow < matrix1.rows(); iRow++)
  {
    for (int iCol = 0; iCol < matrix1.cols(); iCol++)
    {
      if (matrix1(iRow,iCol) < min and matrix2(iRow,iCol) < min) 
      {
        matrix1(iRow,iCol) = 1.0;
        matrix2(iRow,iCol) = 1.0;
      }
    }
  }
  residualMatrix = (ones-(matrix2.cwiseQuotient(matrix1)));
  residual = (1.0/residualMatrix.size())*residualMatrix.squaredNorm();
  return residual;

};
//==============================================================================


//==============================================================================
/// Update the transport fluxes with those in the quasidiffusion objects
void TransportToQDCoupling::updateTransportFluxes()
{
  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    MGT->SGTs[iGroup]->sFlux = MGQD->SGQDs[iGroup]->sFlux;
  }
}
//=============================================================================

//==============================================================================
/// Update transport fluxes at the previous time step with the fluxes currently
/// on the quasidiffusion objects
void TransportToQDCoupling::updateTransportPrevFluxes()
{
  for (int iGroup = 0; iGroup < MGT->SGTs.size(); iGroup++)
  {
    MGT->SGTs[iGroup]->sFluxPrev = MGQD->SGQDs[iGroup]->sFlux;
  }
}
//=============================================================================

//==============================================================================
/// Check for optional input parameters of relevance to this object
void TransportToQDCoupling::checkOptionalParams()
{
  if ((*input)["parameters"]["epsEddington"])
  {
    epsEddington=(*input)["parameters"]["epsEddington"].as<double>();
  }
}
//==============================================================================
