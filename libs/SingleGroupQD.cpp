// File: SingleGroupQD.cpp
// Purpose: define a single group quasidiffusion object
// Date: February 5, 2020

#include "SingleGroupQD.h"

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
SingleGroupQD::SingleGroupQD(int myEnergyGroup,\
    MultiGroupQD * myMGQD,\
    Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * myInput)
{

  // Assign energy group for this object 
  energyGroup = myEnergyGroup;

  // Assign pointers
  MGQD = myMGQD;
  mats = myMaterials;
  mesh = myMesh;
  input = myInput;

  // initialize Eddington factors to diffusion physics
  double diagValue = 1.0/3.0, offDiagValue = 0.0;
  
  Err.setConstant(mesh->nZ,mesh->nR,diagValue);
  Ezz.setConstant(mesh->nZ,mesh->nR,diagValue);
  Erz.setConstant(mesh->nZ,mesh->nR,offDiagValue);
  ErrAxial.setConstant(mesh->nZ+1,mesh->nR,diagValue);
  EzzAxial.setConstant(mesh->nZ+1,mesh->nR,diagValue);
  ErzAxial.setConstant(mesh->nZ+1,mesh->nR,offDiagValue);
  ErrRadial.setConstant(mesh->nZ,mesh->nR+1,diagValue);
  EzzRadial.setConstant(mesh->nZ,mesh->nR+1,diagValue);
  ErzRadial.setConstant(mesh->nZ,mesh->nR+1,offDiagValue);
  G.setConstant(mesh->nZ,mesh->nR,offDiagValue);
  GRadial.setConstant(mesh->nZ,mesh->nR+1,offDiagValue);

  // initialize previous Eddington factors
  ErrPrev = Err;
  EzzPrev = Ezz;
  ErzPrev = Erz;

  // initialize source 
  q.setOnes(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  q = 0.0*q;

  // initialize flux and current matrices
  sFlux.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  sFluxR.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size()+1);
  sFluxZ.setZero(mesh->zCornerCent.size()+1,mesh->rCornerCent.size());
  currentR.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size()+1);
  currentZ.setZero(mesh->zCornerCent.size()+1,mesh->rCornerCent.size());

  // Initialize past values
  sFluxPrev.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  sFluxRPrev.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size()+1);
  sFluxZPrev.setZero(mesh->zCornerCent.size()+1,mesh->rCornerCent.size());
  currentRPrev.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size()+1);
  currentZPrev.setZero(mesh->zCornerCent.size()+1,mesh->rCornerCent.size());

  // initialize boundary conditions
  wFluxBC.setZero(mesh->dzsCorner.size());
  eFluxBC.setZero(mesh->dzsCorner.size());
  nFluxBC.setZero(mesh->drsCorner.size());
  sFluxBC.setZero(mesh->drsCorner.size());
  wCurrentRBC.setZero(mesh->dzsCorner.size());
  eCurrentRBC.setZero(mesh->dzsCorner.size());
  nCurrentZBC.setZero(mesh->drsCorner.size());
  sCurrentZBC.setZero(mesh->drsCorner.size());
  eInwardFluxBC.setZero(mesh->dzsCorner.size());
  nInwardFluxBC.setZero(mesh->drsCorner.size());
  sInwardFluxBC.setZero(mesh->drsCorner.size());
  eInwardCurrentBC.setZero(mesh->dzsCorner.size());
  nInwardCurrentBC.setZero(mesh->drsCorner.size());
  sInwardCurrentBC.setZero(mesh->drsCorner.size());
  eOutwardCurrToFluxRatioBC.setZero(mesh->dzsCorner.size());
  nOutwardCurrToFluxRatioBC.setZero(mesh->drsCorner.size());
  sOutwardCurrToFluxRatioBC.setZero(mesh->drsCorner.size());
  eAbsCurrentBC.setZero(mesh->dzsCorner.size());
  nAbsCurrentBC.setZero(mesh->drsCorner.size());
  sAbsCurrentBC.setZero(mesh->drsCorner.size());

  // check for optional parameters
  checkOptionalParams();
};
//==============================================================================

//==============================================================================
/// Build the portion of low order QD linear system belonging to this object
void SingleGroupQD::formContributionToLinearSystem()
{
  MGQD->QDSolve->formLinearSystem(this);
}
//==============================================================================

//==============================================================================
/// Build the portion of low order QD linear system belonging to this object
void SingleGroupQD::formSteadyStateContributionToLinearSystem()
{
  MGQD->QDSolve->formSteadyStateLinearSystem(this);
}
//==============================================================================

//==============================================================================
/// Build the portion of linear system, used to calculate currents, belonging 
/// to this object
void SingleGroupQD::formContributionToBackCalcSystem()
{
  MGQD->QDSolve->formBackCalcSystem(this);
}
//==============================================================================

//==============================================================================
/// Build the portion of linear system, used to calculate currents, belonging 
/// to this object
void SingleGroupQD::formSteadyStateContributionToBackCalcSystem()
{
  MGQD->QDSolve->formSteadyStateBackCalcSystem(this);
}
//==============================================================================


//==============================================================================
/// Get the flux for this object
void SingleGroupQD::getFlux()
{
  MGQD->QDSolve->getFlux(this);
}
//==============================================================================

//==============================================================================
/// Map the fluxes of this object into a 1D solution vector
Eigen::VectorXd SingleGroupQD::getFluxSolutionVector()
{
  return MGQD->QDSolve->getFluxSolutionVector(this);
}
//==============================================================================

//==============================================================================
/// Map the current of this object into a 1D solution vector
Eigen::VectorXd SingleGroupQD::getCurrentSolutionVector()
{
  return MGQD->QDSolve->getCurrentSolutionVector(this);
}
//==============================================================================

//==============================================================================
/// Print BC parameters 
///
void SingleGroupQD::printEddingtons()
{

cout << "Err: " << endl;
cout << Err << endl;
cout << endl;

cout << "Ezz: " << endl;
cout << Ezz << endl;
cout << endl;

cout << "Erz: " << endl;
cout << Erz << endl;
cout << endl;

cout << "ErrAxial: " << endl;
cout << ErrAxial << endl;
cout << endl;

cout << "EzzAxial: " << endl;
cout << EzzAxial << endl;
cout << endl;

cout << "ErzAxial: " << endl;
cout << ErzAxial << endl;
cout << endl;

cout << "ErrRadial: " << endl;
cout << ErrRadial << endl;
cout << endl;

cout << "EzzRadial: " << endl;
cout << EzzRadial << endl;
cout << endl;

cout << "ErzRadial: " << endl;
cout << ErzRadial << endl;
cout << endl;

};
//==============================================================================


//==============================================================================
// Check for optional input parameters relevant to this object
void SingleGroupQD::checkOptionalParams()
{
  // vectors for reading in optional parameters
  vector<double> inpaFlux,inpCurrent,inpSFluxPrev,inpCurrentPrev;

  // check for optional parameters specified in input file

  if ((*input)["parameters"]["lowerBC"]){

    inpaFlux=(*input)["parameters"]["lowerBC"]\
             .as<vector<double> >();

    nFluxBC.setOnes();
    if (inpaFlux.size() == 1)
      nFluxBC = 4.0*inpaFlux[0]*nFluxBC;
    else
      nFluxBC = 4.0*inpaFlux[energyGroup]*nFluxBC;
  }

  if ((*input)["parameters"]["upperBC"]){

    inpaFlux=(*input)["parameters"]["upperBC"]\
             .as<vector<double> >();

    sFluxBC.setOnes();
    if (inpaFlux.size() == 1)
      sFluxBC = 4.0*inpaFlux[0]*sFluxBC;
    else
      sFluxBC = 4.0*inpaFlux[energyGroup]*sFluxBC;
  }

  if ((*input)["parameters"]["innerBC"]){

    inpaFlux=(*input)["parameters"]["innerBC"]\
             .as<vector<double> >();

    wFluxBC.setOnes();
    if (inpaFlux.size() == 1)
      wFluxBC = 4.0*inpaFlux[0]*wFluxBC;
    else
      wFluxBC = 4.0*inpaFlux[energyGroup]*wFluxBC;
  }

  if ((*input)["parameters"]["outerBC"]){

    inpaFlux=(*input)["parameters"]["outerBC"]\
             .as<vector<double> >();

    eFluxBC.setOnes();
    if (inpaFlux.size() == 1)
      eFluxBC = 4.0*inpaFlux[0]*eFluxBC;
    else
      eFluxBC = 4.0*inpaFlux[energyGroup]*eFluxBC;
  }

  if ((*input)["parameters"]["lowerCurrentBC"]){

    inpCurrent=(*input)["parameters"]["lowerCurrentBC"]\
               .as<vector<double> >();

    nCurrentZBC.setOnes();
    if (inpaFlux.size() == 1)
      nCurrentZBC = inpCurrent[0]*nCurrentZBC;
    else
      nCurrentZBC = inpCurrent[energyGroup]*nCurrentZBC;
  }

  if ((*input)["parameters"]["upperCurrentBC"]){

    inpCurrent=(*input)["parameters"]["upperCurrentBC"]\
               .as<vector<double> >();

    sCurrentZBC.setOnes();
    if (inpaFlux.size() == 1)
      sCurrentZBC = inpCurrent[0]*sCurrentZBC;
    else
      sCurrentZBC = inpCurrent[energyGroup]*sCurrentZBC;
  }

  if ((*input)["parameters"]["innerCurrentBC"]){

    inpCurrent=(*input)["parameters"]["innerCurrentBC"]\
               .as<vector<double> >();

    wCurrentRBC.setOnes();
    if (inpaFlux.size() == 1)
      wCurrentRBC = inpCurrent[0]*wCurrentRBC;
    else
      wCurrentRBC = inpCurrent[energyGroup]*wCurrentRBC;
  }

  if ((*input)["parameters"]["outerCurrentBC"]){

    inpCurrent=(*input)["parameters"]["outerCurrentBC"]\
               .as<vector<double> >();

    eCurrentRBC.setOnes();
    if (inpaFlux.size() == 1)
      eCurrentRBC = inpCurrent[0]*eCurrentRBC;
    else
      eCurrentRBC = inpCurrent[energyGroup]*eCurrentRBC;
  }

  if ((*input)["parameters"]["initial previous flux"]){

    inpSFluxPrev=(*input)["parameters"]["initial previous flux"]\
                 .as<vector<double> >();

    sFlux.setOnes();
    sFluxR.setOnes();
    sFluxZ.setOnes();
    if (inpSFluxPrev.size() == 1)
    {
      sFlux = inpSFluxPrev[0]*sFlux;
      sFluxR = inpSFluxPrev[0]*sFluxR;
      sFluxZ = inpSFluxPrev[0]*sFluxZ;
    }
    else
    {
      sFlux = inpSFluxPrev[energyGroup]*sFlux;
      sFluxR = inpSFluxPrev[energyGroup]*sFluxR;
      sFluxZ = inpSFluxPrev[energyGroup]*sFluxZ;
    }
  }

  if ((*input)["parameters"]["initial previous current"]){

    inpCurrentPrev=(*input)["parameters"]["initial previous current"]\
                   .as<vector<double> >();

    currentR.setOnes();
    currentZ.setOnes();
    if (inpSFluxPrev.size() == 1)
    {
      currentR = inpCurrentPrev[0]*currentR;
      currentZ = inpCurrentPrev[0]*currentZ;
    }
    else
    {
      currentR = inpCurrentPrev[energyGroup]*currentR;
      currentZ = inpCurrentPrev[energyGroup]*currentZ;
    }
  }

}
//==============================================================================

//==============================================================================
/// Write the flux in this energy group to a CVS

void SingleGroupQD::writeFlux()
{

  ofstream fluxFile;
  string fileName;

  fileName = "qd-scalar-flux-group-" + to_string(energyGroup)\
              +"-dz-" + to_string(mesh->dz)+ "-dt-"+to_string(mesh->dt)\
              + "-T-" + to_string(mesh->T) + ".csv";

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

