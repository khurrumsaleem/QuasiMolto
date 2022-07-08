// File: SingleGroupTransport.cpp
// Purpose: define a single group transport object
// Date: October 28, 2019

#include "SingleGroupTransport.h"
#include "MultiGroupTransport.h"
#include "StartingAngle.h"
#include "SimpleCornerBalance.h"

using namespace std; 

//==============================================================================
/// Class object constructor
///
/// @param [in] myEnergyGroup Energy group for this single group transport
/// object
/// @param [in] myMGT MultiGroupTransport object for this simulation
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] myInput YAML input object for the simulation
SingleGroupTransport::SingleGroupTransport(int myEnergyGroup,\
    MultiGroupTransport * myMGT,\
    Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput)
{
  // Assign energy group for this object 
  energyGroup = myEnergyGroup;

  // Assign pointers
  MGT = myMGT;
  mats = myMaterials;
  mesh = myMesh;
  input = myInput;

  // vectors for reading in temporary variables
  vector<double> inpSFlux0,inpSFluxPrev0,inpAlpha0;

  // Initialize angular fluxes
  aFlux.set_size(mesh->zCornerCent.size(),mesh->rCornerCent.size(),mesh->nAngles);
  aFlux.zeros();

  // Initialize half angle angular fluxes used to approximate the angular
  // redistribution term of the RZ neutron transport equation
  aHalfFlux.set_size(mesh->zCornerCent.size(),mesh->rCornerCent.size(),\
      mesh->quadrature.size());
  aHalfFlux.zeros();

  // Initialize scalar fluxes
  sFlux.setOnes(mesh->zCornerCent.size(),mesh->rCornerCent.size());

  // Set scalar flux at previous time step. Set to ones. Setting to zero
  // would produce NaNs when calculating the alphas in each cell
  sFluxPrev.setOnes(mesh->zCornerCent.size(),mesh->rCornerCent.size());

  // Initialize alphas, total source, scattering source, and fission source
  alpha.setOnes(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  q.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  scatterSource.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  fissionSource.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size());

  // Check for optional inputs and read them in.
  if ((*input)["parameters"]["initial flux"]){

    inpSFlux0=(*input)["parameters"]["initial flux"].as<vector<double>>();

    if (inpSFlux0.size() == 1)
      sFlux = inpSFlux0[0]*sFlux;
    else
      sFlux = inpSFlux0[energyGroup]*sFlux;
  } 

  if ((*input)["parameters"]["initial previous flux"]){

    inpSFluxPrev0=(*input)["parameters"]["initial previous flux"].as<vector<double>>();

    if (inpSFluxPrev0.size() == 1)
      sFluxPrev = inpSFluxPrev0[0]*sFluxPrev;
    else
      sFluxPrev = inpSFluxPrev0[energyGroup]*sFluxPrev;
  } 

  if ((*input)["parameters"]["initial alpha"]){

    inpAlpha0=(*input)["parameters"]["initial alpha"].as<vector<double>>();

    if (inpAlpha0.size() == 1)
      alpha = inpAlpha0[0]*alpha;
    else
      alpha = inpAlpha0[energyGroup]*alpha;
  } 

};

//==============================================================================

//==============================================================================
/// Call starting angle solver

void SingleGroupTransport::solveStartAngle()
{
  aHalfFlux.zeros();
  MGT->startAngleSolve->calcStartingAngle(&aHalfFlux,&q,&alpha,energyGroup);
};

//==============================================================================

//==============================================================================
/// Call simple corner balance solver

void SingleGroupTransport::solveSCB()
{
  aFlux.zeros();  
  MGT->SCBSolve->solve(&aFlux,&aHalfFlux,&q,&alpha,energyGroup);
};

//==============================================================================

//==============================================================================
/// Calculate the source in this object
///
/// Assuming the fluxes currently contained in this object
/// @param [in] calcType A string that indicates how to compute the source. 
/// "s" indicates to recalculate the scatter source, and not recalculate the
/// fission source. "fs" (default) recalculate the fission and scattering 
/// source.
/// @param [out] residual L2 norm of difference between the newly calculated 
/// and past source
double SingleGroupTransport::calcSource(string calcType)
{ 
  double residual; 
  Eigen::MatrixXd q_old = q;

  // Check is grey group sources should be used
  if (MGT->useMPQDSources)
  {         
    if (calcType == "steady_state")
      q = calcSteadyStateMPQDSource();
    else
      q = calcMPQDSource();
  }
  // Otherwise, use transport sources
  else
  {
    if (calcType=="s" or calcType=="S")
    {
      // Just re-evaluate scattering source
      calcScatterSource();
      q = scatterSource + fissionSource;
    } 
    else if (calcType == "fs" or calcType == "FS")
    {
      // Re-evaluate scattering and fission source
      calcScatterSource();
      calcFissionSource();
      q = scatterSource + fissionSource;
    }
  }

  // Calculate residual
  residual =  (q_old-q).norm();

  return residual;
};

//==============================================================================

//==============================================================================
/// Calculate the fission source in this object
///
/// Assuming the fluxes currently contained in this object
/// @param [out] residual L2 norm of difference between the newly calculated 
/// fission source and the past fission source
double SingleGroupTransport::calcFissionSource()
{  

  double residual;
  Eigen::MatrixXd fissionSource_old = fissionSource;
  Eigen::MatrixXd fissionSourceDiff;

  // Set fission source to zero
  fissionSource.setZero();

  // Calculate fission source 
  for (int iZ = 0; iZ < fissionSource.rows(); ++iZ){
    for (int iR = 0; iR < fissionSource.cols(); ++iR){
      for (int iGroup = 0; iGroup < MGT->SGTs.size(); ++iGroup){

        fissionSource(iZ,iR) = fissionSource(iZ,iR) + (1.0/mesh->totalWeight)*(
            mats->chiP(iZ,iR,energyGroup)*mats->nu(iZ,iR,energyGroup)\
            *mats->sigF(iZ,iR,iGroup)*MGT->SGTs[iGroup]->sFlux(iZ,iR));

        // need to account for precursors, too.
      } // iGroup
    } // iR
  } // iZ

  fissionSourceDiff = (fissionSource_old-fissionSource);

  // Need to be careful how we calculate the relative difference here, as
  // some values we're dividing by may be equal to 0.
  for (int iRow = 0; iRow < fissionSourceDiff.rows(); ++iRow){
    for (int iCol = 0; iCol < fissionSourceDiff.cols(); ++iCol){
      if (fissionSource(iRow,iCol)!=0){
        fissionSourceDiff(iRow,iCol) = fissionSourceDiff(iRow,iCol)\
                                       /fissionSource(iRow,iCol);
      } else {
        fissionSourceDiff(iRow,iCol) = 0;
      }
    }
  }

  // Calculate residual
  residual=fissionSourceDiff.norm();
  return residual;
};

//==============================================================================

//==============================================================================
/// Calculate the scatter source in this object
///
/// Assuming the fluxes currently contained in this object
/// @param [out] residual L2 norm of difference between newly calculated 
/// scattering source and the past scattering source
double SingleGroupTransport::calcScatterSource()
{  
  double residual;
  Eigen::MatrixXd scatterSource_old = scatterSource;

  // Set scattering source to zero
  scatterSource.setZero();

  // Calculate the scattering source
  for (int iZ = 0; iZ < scatterSource.rows(); ++iZ){
    for (int iR = 0; iR < scatterSource.cols(); ++iR){
      for (int iGroup = 0; iGroup < MGT->SGTs.size(); ++iGroup){
        scatterSource(iZ,iR) = scatterSource(iZ,iR) + (1.0/mesh->totalWeight)*(
            mats->sigS(iZ,iR,iGroup,energyGroup)*MGT->SGTs[iGroup]->sFlux(iZ,iR));

        // need to account for precursors, too.
      } // iGroup
    } // iR
  } // iZ

  // Calculate residual
  residual = ((scatterSource_old-scatterSource)\
      .cwiseQuotient(scatterSource)).norm();

  return residual;
};

//==============================================================================

//==============================================================================
/// Calculate the grey group source 
///
Eigen::MatrixXd SingleGroupTransport::calcMPQDSource()
{  
  double localSigS, localChiP, localChiD, localFissionCoeff, localFlux,\
    localDNPSource;
  double weight = 1.0/mesh->totalWeight; 
  Eigen::MatrixXd mpqdSource;

  // Initialize size of source matrix
  mpqdSource.setZero(sFlux.rows(),sFlux.cols());

  for (int iR = 0; iR < mpqdSource.cols(); iR++)
  {
    for (int iZ = 0; iZ < mpqdSource.rows(); iZ++)
    {
      localSigS = mats->oneGroupXS->groupScatterXS(iZ,iR,energyGroup);
      localFissionCoeff = mats->oneGroupXS->qdFluxCoeff(iZ,iR);
      localChiP = mats->chiP(iZ,iR,energyGroup);
      localChiD = mats->chiD(iZ,iR,energyGroup);
      localFlux = MGT->mpqd->ggqd->sFlux(iZ,iR); 
      localDNPSource = MGT->mpqd->mgdnp->dnpSource(iZ,iR); 

      // Scattering source
      mpqdSource(iZ,iR) = weight*localSigS*localFlux;

      // Fission source
      mpqdSource(iZ,iR) += weight*localChiP*localFissionCoeff*localFlux;

      // DNP source
      mpqdSource(iZ,iR) += weight*localChiD*localDNPSource;
    }
  }

  return mpqdSource;
};

//==============================================================================

//==============================================================================
/// Calculate the steady state grey group source 
///
Eigen::MatrixXd SingleGroupTransport::calcSteadyStateMPQDSource()
{  
  double localSigS, localChiP, localChiD, localFissionCoeff, localFlux,\
    localDNPSource, keff;
  double weight = 1.0 / mesh->totalWeight; 
  Eigen::MatrixXd mpqdSource;

  // Initialize size of source matrix
  mpqdSource.setZero(sFlux.rows(),sFlux.cols());

  for (int iR = 0; iR < mpqdSource.cols(); iR++)
  {
    for (int iZ = 0; iZ < mpqdSource.rows(); iZ++)
    {
      localSigS = mats->oneGroupXS->groupScatterXS(iZ,iR,energyGroup);
      localFissionCoeff = mats->oneGroupXS->qdFluxCoeff(iZ,iR);
      localChiP = mats->chiP(iZ,iR,energyGroup);
      localChiD = mats->chiD(iZ,iR,energyGroup);
      localFlux = MGT->mpqd->ggqd->sFlux(iZ,iR); 
      localDNPSource = MGT->mpqd->mgdnp->dnpSource(iZ,iR); 
      keff = mats->oneGroupXS->keff;

      // Scattering source
      mpqdSource(iZ,iR) = weight*localSigS*localFlux;

      // Fission source
      mpqdSource(iZ,iR) += weight*localChiP*localFissionCoeff*localFlux/keff;

      // DNP source
      mpqdSource(iZ,iR) += weight*localChiD*localDNPSource;
    }
  }

  return mpqdSource;
};

//==============================================================================

//==============================================================================
/// Calculate the flux in this energy group
///
/// Assuming the fluxes currently contained in this object
/// @param [out] residual L2 norm of the difference between the newly calculated
/// flux and past flux
double SingleGroupTransport::calcFlux()
{
  double weight,residual; 
  int weightIdx=3,angIdx;
  Eigen::MatrixXd sFlux_old = sFlux;

  // Set scalar flux to zero
  sFlux.setZero();

  // Calculate scalar flux
  for (int iQ = 0; iQ < mesh->quadrature.size(); ++iQ){
    for (int iP = 0; iP < mesh->quadrature[iQ].nOrd; ++iP){

      weight = mesh->quadrature[iQ].quad[iP][weightIdx];
      angIdx = mesh->quadrature[iQ].ordIdx[iP];

      for (int iZ = 0; iZ < sFlux.rows(); ++iZ){
        for (int iR = 0; iR < sFlux.cols(); ++iR){

          sFlux(iZ,iR) = sFlux(iZ,iR)\
                         +weight*aFlux(iZ,iR,angIdx);

        } // iR
      } // iZ
    } // iP
  } // iQ

  // Calculate residual
  residual = ((sFlux_old-sFlux).cwiseQuotient(sFlux)).norm();
  
  return residual;
};

//==============================================================================

//==============================================================================
/// Calculate alpha, which approximates time dependent behaviour, in each cell 
///
/// Assuming the fluxes currently contained in this object
/// @param [out] residual L2 norm of difference between newly calculate alpha
/// and past alpha
double SingleGroupTransport::calcAlpha(string calcType)
{  
  double localFlux,localFluxPrev,residual,deltaT = mesh->dt;
  vector<int> indices;
  Eigen::MatrixXd alpha_old = alpha;
  Eigen::MatrixXd alphaDiff;
  PetscErrorCode ierr = 0;
  PetscScalar value;
  PetscInt index;
  Vec            x_p_seq;
  VecScatter     ctx;

  // Set alpha to zero
  alpha.setZero();

  if (MGT->useMPQDSources)
  {
    if (calcType == "steady_state")
      alpha.setZero();
    else
    {
      VecScatterCreateToAll(MGT->mgqd->QDSolve->x_p,&ctx,&x_p_seq);
      VecScatterBegin(ctx,MGT->mgqd->QDSolve->x_p,x_p_seq,INSERT_VALUES,SCATTER_FORWARD);
      VecScatterEnd(ctx,MGT->mgqd->QDSolve->x_p,x_p_seq,INSERT_VALUES,SCATTER_FORWARD);

      for (int iZ = 0; iZ < alpha.rows(); ++iZ){
        for (int iR = 0; iR < alpha.cols(); ++iR){

          index = MGT->mgqd->QDSolve->getIndices(iR, iZ, energyGroup)[0]; 

          // Get local values
          ierr = VecGetValues(x_p_seq, 1, &index, &localFlux); CHKERRQ(ierr);
          ierr = VecGetValues(MGT->mgqd->QDSolve->xPast_p_seq, 1, &index, &localFluxPrev); CHKERRQ(ierr);
          alpha(iZ, iR) = (1.0 / deltaT) * log(localFlux / localFluxPrev);

        } // iR
      } // iZ

      VecScatterDestroy(&ctx);
      VecDestroy(&x_p_seq);
    }
  }
  else
  {
    // Calculate alphas
    for (int iZ = 0; iZ < alpha.rows(); ++iZ)
    {
      for (int iR = 0; iR < alpha.cols(); ++iR)
      {
        alpha(iZ, iR) = (1.0 / deltaT) * log(sFlux(iZ, iR) / sFluxPrev(iZ, iR));
      } // iR
    } // iZ
  }

  alphaDiff = (alpha_old-alpha);
  // Calculate residual
  // Need to be careful how we calculate the relative difference here, as
  // some values we're dividing by may be equal to 0.
  for (int iRow = 0; iRow < alphaDiff.rows(); ++iRow){
    for (int iCol = 0; iCol < alphaDiff.cols(); ++iCol){
      if (alpha(iRow,iCol)!=0){
        alphaDiff(iRow,iCol) = alphaDiff(iRow,iCol)/alpha(iRow,iCol);
      } else {
        alphaDiff(iRow,iCol) = 0;
      }
    }
  }

  residual = alphaDiff.norm();

  return residual;
};
//==============================================================================

//==============================================================================
/// Write the flux in this energy group to a CVS

void SingleGroupTransport::writeFlux()
{
  ofstream fluxFile;
  string fileName;

  // parse file name 
  fileName = "scalar-flux-group-" + to_string(energyGroup)\
              +"-dz-" + to_string(mesh->dz)+ "-dt-"+to_string(mesh->dt)\
              +"-T-"+ to_string(mesh->T) + ".csv";

  // open file
  fluxFile.open(fileName);

  // write flux values to .csv
  for (int iZ = 0; iZ < sFlux.rows(); ++iZ) {
    fluxFile << sFlux(iZ,0);
    for (int iR = 1; iR < sFlux.cols(); ++iR) {
      fluxFile <<","<< sFlux(iZ,iR);
    }
    fluxFile << endl;
  }
  fluxFile.close();

  // parse file name 
  fileName = "angular-flux-group-" + to_string(energyGroup)\
              +"-dz-" + to_string(mesh->dz)+ "-dt-"+to_string(mesh->dt)\
              +"-T-"+ to_string(mesh->T) + ".csv";

  // open file
  fluxFile.open(fileName);

  // write flux values to .csv
  for (int iZ = 0; iZ < sFlux.rows(); ++iZ) {
    fluxFile << aFlux(iZ,0,0);
    for (int iR = 1; iR < sFlux.cols(); ++iR) {
      fluxFile <<","<< aFlux(iZ,iR,0);
    }
    fluxFile << endl;
  }
  fluxFile.close();

  // if this is the first energy group, write mesh too
  if (energyGroup==0){
    // write radial mesh to .csv 
    fileName = "r-mesh.csv"; 
    fluxFile.open(fileName); 
    fluxFile << mesh->rCornerCent(0);
    for (int iR = 1; iR < mesh->rCornerCent.size(); ++iR) {
      fluxFile << "," << mesh->rCornerCent(iR);
    }
    fluxFile << endl;
    fluxFile.close();

    // write axial mesh to .csv
    fileName = "z-mesh.csv"; 
    fluxFile.open(fileName); 
    fluxFile << mesh->zCornerCent(0);
    for (int iZ = 1; iZ < mesh->zCornerCent.size(); ++iZ) {
      fluxFile << "," << mesh->zCornerCent(iZ);
    }
    fluxFile << endl;
    fluxFile.close();
  }
};
//==============================================================================
