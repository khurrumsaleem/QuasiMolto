// File: Materials.cpp
// Purpose: contain material characteristics and geometry for a simulation
// Date: October 28, 2019

#include "Materials.h"

using namespace std; 

//==============================================================================
/// Material class object constructor
///
/// @param [in] myMesh Mesh object for this simulation
/// @param [in] myInput YAML input file for this simulation
Materials::Materials(Mesh * myMesh,YAML::Node * myInput)
{
  // Temporary vector for holding neutron velocities before reading them into
  // and Eigen data type
  vector<double> neutVInp;    

  // Point to variables for mesh and input file
  mesh = myMesh;
  input = myInput;

  // Initialize material map
  matMap.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size());

  // Initialize flow velocity matrix in core 
  flowVelocity.setZero(mesh->zCornerCent.size(),mesh->rCornerCent.size());
  
  // Initialize temperature in core 
  temperature.setConstant(mesh->nZ,mesh->nR,922.0);

  // Initialize flow velocity matrix in recirculation loop
  recircFlowVelocity.setZero(mesh->nZrecirc,mesh->rCornerCent.size());

  // Read materials from input
  readMats();

  // Read geometry information from input
  readGeom();

  // Read flow velocity information from input
  readFlowVelocity();

  // Print information on materials and geometry, and check to see if nuclear 
  // data is sensible
  edit();
  checkMats();        

  // Check if neutron velocities are specified in input
  if ((*input)["parameters"]["neutron velocity"]){


    neutVInp = (*input)["parameters"]["neutron velocity"]\
               .as<vector<double>>();

    if (neutVInp.size() == 1){
      neutV.setOnes(nGroups);
      neutV = neutVInp[0]*neutV;
    } 
    else {

      neutV.setZero(nGroups);
      for (int iGroup = 0; iGroup < nGroups; ++iGroup){
        neutV(iGroup) = neutVInp[iGroup];
      }
    }
  } else {

    // Set default neutron velocities 
    neutV.setOnes(nGroups);
    neutV = 2200.0*neutV;
  }

  if ((*input)["parameters"]["uniformTemperature"]){
    uniformTemperature = (*input)["parameters"]["uniformTemperature"]\
               .as<bool>();
    uniformTempValue = (*input)["parameters"]["inletTemp"]\
               .as<double>();
  }

  // Initialize 1GXS
  oneGroupXS = new CollapsedCrossSections(mesh,nGroups);

  // Initialize data in collapsed cross sections class
  initCollapsedXS();

};
//==============================================================================

//==============================================================================
/// Parse input and define geometry parameters

void Materials::readGeom()
{
  YAML::Node geometry = (*input)["geometry"];
  string regType,matName;
  double rIn,rOut,zLow,zUp;
  int index;

  // Loop over material section  
  for (YAML::const_iterator it=geometry.begin(); it!=geometry.end(); ++it){

    // Check material type
    regType=it->first.as<string>();

    if (regType=="background"){

      // Set the entire material map to this material
      index = mat2idx[it->second.as<string>()];
      setMatRegion(index,0.0,mesh->R,0.0,mesh->Z);

    } else if (regType=="region"){

      // Set a portion of material map to this material
      index = mat2idx[it->second["material"].as<string>()];

      // Get bounds of this material
      rIn = it->second["inner-r"].as<double>(); 
      rOut = it->second["outer-r"].as<double>(); 
      zLow = it->second["lower-z"].as<double>(); 
      zUp = it->second["upper-z"].as<double>();

      // Make call set region index
      setMatRegion(index,rIn,rOut,zLow,zUp);
    }
  }
};
//==============================================================================

//==============================================================================
/// Parse input and define material parameters

void Materials::readMats()
{
  YAML::Node mats = (*input)["materials"];
  string name;
  vector<double> sigTInp,sigSInp,sigFInp,chiPInp,chiDInp,nuInp,neutVInp;
  vector<Eigen::MatrixXd> sigT,sigS,sigF,nu,neutV;
  Eigen::VectorXd chiP,chiD; 
  int ID,size;
  double density,gamma,k,cP,omega,flowVelocity;
  bool stationary;

  int iCount=0;
  for (YAML::const_iterator it=mats.begin(); it!=mats.end(); ++it){
    cout << iCount << endl;
    
    // reset boolean
    stationary = true;

    // Add material to material to index map  
    mat2idx.insert(pair<string,int>(it->first.as<string>(),iCount));
    name = it->first.as<string>();

    // Pull sigT values from input
    if (it->second["sigTFile"])
      sigT=readTempDependentYaml(it->second["sigTFile"].as<string>());
    else
    {
      sigTInp = it->second["sigT"].as<vector<double>>();
      sigT.resize(sigTInp.size()); 

      for (int iSig = 0; iSig < sigT.size(); ++iSig){
        sigT[iSig].setZero(1,2);
        sigT[iSig](0,1) = sigTInp[iSig];
      }
    }
    size = sigT.size();

    // Pull sigF values from input
    if (it->second["sigFFile"])
      sigF=readTempDependentYaml(it->second["sigFFile"].as<string>());
    else
    {
      sigFInp = it->second["sigF"].as<vector<double>>();
      sigF.resize(size); 

      for (int iSig = 0; iSig < size; ++iSig){
        sigF[iSig].setZero(1,2);
        sigF[iSig](0,1) = sigFInp[iSig];
      }
    }

    // Pull sigS values from input
    if (it->second["sigSFile"])
      sigS=readTempDependentYaml(it->second["sigSFile"].as<string>());
    else
    {
      sigSInp = it->second["sigS"].as<vector<double>>();
      sigS.resize(sigSInp.size());

      for (int iSig = 0; iSig < size; ++iSig){
        
        for(int iGroup = 0; iGroup < size; ++iGroup){
          sigS[iSig*size+iGroup].setZero(1,2);
          sigS[iSig*size+iGroup](0,1) = sigSInp[iSig*size+iGroup];
        }
      }
    }
    
    // Pull nu values from input
    if (it->second["nuFile"])
      nu = readTempDependentYaml(it->second["nuFile"].as<string>());
    else
    {
      nuInp = it->second["nu"].as<vector<double>>();
      nu.resize(size); 

      for (int iNu = 0; iNu < size; ++iNu){
        nu[iNu].setZero(1,2);
        nu[iNu](0,1) = nuInp[iNu];
      }
    }

    // Pull neutron velocity values from input
    if (it->second["neutVFile"])
      neutV = readTempDependentYaml(it->second["neutVFile"].as<string>());
    else
    {
      neutVInp = it->second["neutV"].as<vector<double>>();
      neutV.resize(size); 

      for (int iNeut = 0; iNeut < size; ++iNeut){
        neutV[iNeut].setZero(1,2);
        neutV[iNeut](0,1) = neutVInp[iNeut];
      }
    }


    chiPInp = it->second["chiP"].as<vector<double>>();
    chiDInp = it->second["chiD"].as<vector<double>>();
    density = it->second["density"].as<double>();
    gamma = it->second["gamma"].as<double>();
    k = it->second["k"].as<double>();
    cP = it->second["cP"].as<double>();
    omega = it->second["omega"].as<double>();
    
    if (it->second["material velocity"]) 
      flowVelocity = it->second["material velocity"].as<double>();
    
    if (it->second["stationary"]) 
      stationary = it->second["stationary"].as<bool>();

    // Set size of Eigen vectors
    chiP.setZero(size);
    chiD.setZero(size);

    // Load vector inputs into Eigen vectors
    for (int iSig = 0; iSig < size; ++iSig){

      chiP(iSig) = chiPInp[iSig];
      chiD(iSig) = chiDInp[iSig];

    }

    // Add material to bank
    shared_ptr<Material> newMat (new Material(iCount,name,sigT,sigS,\
          sigF,nu,neutV,chiP,chiD,density,gamma,k,cP,omega,flowVelocity,\
          stationary));
    matBank.push_back(std::move(newMat));
    ++iCount;
  }

  // Set number of energy groups for the simulation	
  nGroups = size;
};
//==============================================================================

//==============================================================================
/// Build matrix describing flow rate in each cell 
///
void Materials::readFlowVelocity()
{

  double inputFlowVelocity;
  Eigen::VectorXd tempVelocity;

  if ((*input)["parameters"]["uniform flow velocity"])
    inputFlowVelocity = (*input)["parameters"]["uniform flow velocity"]\
                        .as<double>();
  else
    inputFlowVelocity = 0.0;

  for (int iZ = 0; iZ < mesh->zCornerCent.size(); iZ++){
    for (int iR = 0; iR < mesh->rCornerCent.size(); iR++){

      // If the material here is not stationary, assign the input flow velocity
      if (!matBank[matMap(iZ,iR)]->stationary){

        flowVelocity(iZ,iR) = coreFlowVelocity(iZ,iR);

      }
    }
  }

  // Check whether velocities are all positive 
  posVelocity = -1E-10 < flowVelocity.minCoeff();

  cout << flowVelocity << endl;
  cout << posVelocity << endl;

  for (int iR = 0; iR < recircFlowVelocity.cols(); ++iR)
  {
    if (posVelocity)
    {
      tempVelocity.setConstant(mesh->nZrecirc,flowVelocity(Eigen::last,iR));
      recircFlowVelocity(Eigen::all,iR) = tempVelocity;
    }
    else
    {
      tempVelocity.setConstant(mesh->nZrecirc,flowVelocity(0,iR));
      recircFlowVelocity(Eigen::all,iR) = tempVelocity;
    }
  }
  cout << "Recirc flow velocity" << endl;
  cout << recircFlowVelocity << endl;
};
//==============================================================================

//==============================================================================
/// Set material region to provided index
///
/// Set all locations within the input bounds equal to the input index
/// @param [in] myIndex Material index to set region to
/// @param [in] rIn Inner radial boundary
/// @param [in] rOut Outer radial boundary
/// @param [in] zLow Lower axial boundary
/// @param [in] zUp Upper axial boundary
void Materials::setMatRegion(int myIndex,double rIn,double rOut,double zLow,double zUp)
{
  for (int iZ = 0; iZ < mesh->zCornerCent.size(); iZ++){
    if (mesh->zCornerCent(iZ)>zLow && mesh->zCornerCent(iZ)<zUp){
      for (int iR = 0; iR < mesh->rCornerCent.size(); iR++){
        if (mesh->rCornerCent(iR)>rIn && mesh->rCornerCent(iR)<rOut){

          matMap(iZ,iR) = myIndex;

        }
      }
    }
  }
};
//==============================================================================

//==============================================================================
/// sigT return total cross section
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] sigT Total cross section of the inquired location 
double Materials::sigT(int zIdx,int rIdx,int eIdx){

  double temp = temperature(zIdx,rIdx);
  double sigT = matBank[matMap(zIdx,rIdx)]->getSigT(eIdx,temp);

  return sigT;
};
//==============================================================================

//==============================================================================
/// sigT return total cross section at axial interface 
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] sigT Total cross section of the inquired location 
double Materials::zSigT(int zIdx,int rIdx,int eIdx){

  double tempDown, tempUp, dzDown, dzUp, sigTDown, sigTUp, sigT;

  if (zIdx == 0)
  {
    tempDown = temperature(zIdx,rIdx);
    sigT = matBank[matMap(zIdx,rIdx)]->getSigT(eIdx,tempDown);
  }
  else if (zIdx == mesh->nZ)
  {
    tempUp = temperature(zIdx-1,rIdx);
    sigT = matBank[matMap(zIdx-1,rIdx)]->getSigT(eIdx,tempUp);
  }
  else
  {
    dzDown = mesh->dzsCorner(zIdx-1);
    tempDown = temperature(zIdx-1,rIdx);
    sigTDown = matBank[matMap(zIdx-1,rIdx)]->getSigT(eIdx,tempDown);
    
    dzUp = mesh->dzsCorner(zIdx);
    tempUp = temperature(zIdx,rIdx);
    sigTUp = matBank[matMap(zIdx,rIdx)]->getSigT(eIdx,tempUp);

    sigT = (sigTDown*dzDown + sigTUp*dzUp)/(dzDown + dzUp);
  }     

  return sigT;
};
//==============================================================================

//==============================================================================
/// sigT return total cross section at radial interface 
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] sigT Total cross section of the inquired location 
double Materials::rSigT(int zIdx,int rIdx,int eIdx){

  double tempLeft, tempRight, drLeft, drRight, sigTLeft, sigTRight, sigT;

  if (rIdx == 0)
  {
    tempLeft = temperature(zIdx,rIdx);
    sigT = matBank[matMap(zIdx,rIdx)]->getSigT(eIdx,tempLeft);
  }
  else if (rIdx == mesh->nR)
  {
    tempRight = temperature(zIdx,rIdx-1);
    sigT = matBank[matMap(zIdx,rIdx-1)]->getSigT(eIdx,tempRight);
  }
  else
  {
    drLeft = mesh->drsCorner(rIdx-1);
    tempLeft = temperature(zIdx,rIdx-1);
    sigTLeft = matBank[matMap(zIdx,rIdx-1)]->getSigT(eIdx,tempLeft);
    
    drRight = mesh->drsCorner(rIdx);
    tempRight = temperature(zIdx,rIdx);
    sigTRight = matBank[matMap(zIdx,rIdx)]->getSigT(eIdx,tempRight);

    sigT = (sigTLeft*drLeft + sigTRight*drRight)/(drLeft + drRight);
  }     

  return sigT;
};
//==============================================================================

//==============================================================================
/// sigS return scattering cross section
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] sigS Scattering cross section of the inquired location 
double Materials::sigS(int zIdx,int rIdx,int gprime,int g){

  double temp = temperature(zIdx,rIdx);
  double sigS = matBank[matMap(zIdx,rIdx)]->getSigS(gprime,g,temp);

  return sigS;
};
//==============================================================================

//==============================================================================
/// sigF return fission cross section 
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] sigF Fission cross section of the inquired location 
double Materials::sigF(int zIdx,int rIdx,int eIndx){

  double temp = temperature(zIdx,rIdx);
  double sigF = matBank[matMap(zIdx,rIdx)]->getSigF(eIndx,temp);

  return sigF;
};
//==============================================================================

//==============================================================================
/// chiP return probability that a prompt neutron will be born 
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] chiP ChiP at the inquired location 
double Materials::chiP(int zIdx,int rIdx,int eIndx){

  double chiP = matBank[matMap(zIdx,rIdx)]->chiP(eIndx);

  // Eventually there will need to be some manipulation here that 
  // extrapolates the cross section based on temperature

  return chiP;
};
//==============================================================================

//==============================================================================
/// chiD return probability that a delayed neutron will be born
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] chiD ChiD at the inquired location 
double Materials::chiD(int zIdx,int rIdx,int eIndx){

  double chiD = matBank[matMap(zIdx,rIdx)]->chiD(eIndx);

  // Eventually there will need to be some manipulation here that 
  // extrapolates the cross section based on temperature

  return chiD;
};
//==============================================================================

//==============================================================================
/// Return nu, the average number of neutrons produced per fission event
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] nu Nu at the inquired location 
double Materials::nu(int zIdx,int rIdx,int eIndx){

  double temp = temperature(zIdx,rIdx);
  double nu = matBank[matMap(zIdx,rIdx)]->getNu(eIndx,temp);

  return nu;
};
//==============================================================================

//==============================================================================
/// Return neutron velocity 
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] neutV Neutron velocity at location and energy
double Materials::neutVel(int zIdx,int rIdx,int eIndx){

  double temp = temperature(zIdx,rIdx);
  double neutV = matBank[matMap(zIdx,rIdx)]->getNeutV(eIndx,temp);

  return neutV;
};
//==============================================================================

//==============================================================================
/// return neutron velocity at axial interface 
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] sigT Total cross section of the inquired location 
double Materials::zNeutVel(int zIdx,int rIdx,int eIdx){

  double tempDown, tempUp, dzDown, dzUp, neutVelDown, neutVelUp, neutVel;

  if (zIdx == 0)
  {
    tempDown = temperature(zIdx,rIdx);
    neutVel = matBank[matMap(zIdx,rIdx)]->getNeutV(eIdx,tempDown);
  }
  else if (zIdx == mesh->nZ)
  {
    tempUp = temperature(zIdx-1,rIdx);
    neutVel = matBank[matMap(zIdx-1,rIdx)]->getNeutV(eIdx,tempUp);
  }
  else
  {
    dzDown = mesh->dzsCorner(zIdx-1);
    tempDown = temperature(zIdx-1,rIdx);
    neutVelDown = matBank[matMap(zIdx-1,rIdx)]->getNeutV(eIdx,tempDown);
    
    dzUp = mesh->dzsCorner(zIdx);
    tempUp = temperature(zIdx,rIdx);
    neutVelUp = matBank[matMap(zIdx,rIdx)]->getNeutV(eIdx,tempUp);

    neutVel = (neutVelDown*dzDown + neutVelUp*dzUp)/(dzDown + dzUp);
  }     

  return neutVel;
};
//==============================================================================

//==============================================================================
/// sigT return total cross section at radial interface 
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] sigT Total cross section of the inquired location 
double Materials::rNeutVel(int zIdx,int rIdx,int eIdx){

  double tempLeft, tempRight, drLeft, drRight, neutVelLeft, neutVelRight,\
    neutVel;

  if (rIdx == 0)
  {
    tempLeft = temperature(zIdx,rIdx);
    neutVel = matBank[matMap(zIdx,rIdx)]->getNeutV(eIdx,tempLeft);
  }
  else if (rIdx == mesh->nR)
  {
    tempRight = temperature(zIdx,rIdx-1);
    neutVel = matBank[matMap(zIdx,rIdx-1)]->getNeutV(eIdx,tempRight);
  }
  else
  {
    drLeft = mesh->drsCorner(rIdx-1);
    tempLeft = temperature(zIdx,rIdx-1);
    neutVelLeft = matBank[matMap(zIdx,rIdx-1)]->getNeutV(eIdx,tempLeft);
    
    drRight = mesh->drsCorner(rIdx);
    tempRight = temperature(zIdx,rIdx);
    neutVelRight = matBank[matMap(zIdx,rIdx)]->getNeutV(eIdx,tempRight);

    neutVel = (neutVelLeft*drLeft + neutVelRight*drRight)/(drLeft + drRight);
  }     

  return neutVel;
};
//==============================================================================

//==============================================================================
/// Return density at indexed locaiton
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] density at inquired location 
double Materials::density(int zIdx,int rIdx){

  double density = matBank[matMap(zIdx,rIdx)]->density;

  return density;
};
//==============================================================================

//==============================================================================
/// Return gamma at indexed locaiton
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] gamma at inquired location 
double Materials::gamma(int zIdx,int rIdx){

  double gamma = matBank[matMap(zIdx,rIdx)]->gamma;

  return gamma;
};
//==============================================================================

//==============================================================================
/// Return k at indexed locaiton
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] k at inquired location 
double Materials::k(int zIdx,int rIdx){

  double k = matBank[matMap(zIdx,rIdx)]->k;

  return k;
};
//==============================================================================

//==============================================================================
/// Return cP at indexed locaiton
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] cP at inquired location 
double Materials::cP(int zIdx,int rIdx){

  double cP = matBank[matMap(zIdx,rIdx)]->cP;

  return cP;
};
//==============================================================================

//==============================================================================
/// Return omega at indexed locaiton
///
/// @param [in] zIdx Z index of location 
/// @param [in] zIdx R index of location 
/// @param [in] eIndx Energy index of location 
/// @param [out] omega at inquired location 
double Materials::omega(int zIdx,int rIdx){

  double omega = matBank[matMap(zIdx,rIdx)]->omega;

  return omega;
};
//==============================================================================

//==============================================================================
/// Return flow velocity at indexed locaiton
///
/// @param [in] zIdx Z index of location 
/// @param [in] rIdx R index of location 
double Materials::coreFlowVelocity(int zIdx,int rIdx){

  double flowVelocity = matBank[matMap(zIdx,rIdx)]->flowVelocity;

  return flowVelocity;
};
//==============================================================================

//==============================================================================
/// Update temperatures used to evaluate cross sections 
///
void Materials::updateTemperature(Eigen::MatrixXd myTemp)
{
  if (not uniformTemperature)
    temperature = myTemp;
  else
    temperature.setConstant(uniformTempValue);
};
//==============================================================================

//==============================================================================
/// Read in temperature dependent data from a YAML file 
///
vector<Eigen::MatrixXd> Materials::readTempDependentYaml(string fileName)
{

  vector<Eigen::MatrixXd> params;
  vector<double> myTemps,myXSs;
  double myTemp,myXS;
  Eigen::MatrixXd myMatrix;

  YAML::Node input;
  input = YAML::LoadFile(fileName);
  for (YAML::const_iterator energy=input.begin();energy!=input.end();++energy) 
  {
    myTemps = energy->second["temperature"].as<vector<double>>();
    myXSs = energy->second["data"].as<vector<double>>();
    myMatrix.setZero(myTemps.size(),2);

    for (int iRow = 0; iRow < myMatrix.rows(); iRow++)
    {
      myMatrix(iRow,0) = myTemps[iRow];
      myMatrix(iRow,1) = myXSs[iRow];
    }
    cout << myMatrix << endl;
    params.push_back(myMatrix);
  }
  return params;
};
//==============================================================================


//==============================================================================
/// Check to see if nuclear data defined on each group is sensible

void Materials::checkMats()
{
  for (int iCount = 0; iCount < matBank.size(); ++iCount){
    matBank[iCount]->checkMat();
  }
};
//==============================================================================

//==============================================================================
/// Initialize collapsed xs data, arbitrarily, to that defined in the fast group 

void Materials::initCollapsedXS()
{

  for (int iR = 0; iR < mesh->nR; ++iR)
  {
    for (int iZ = 0; iZ < mesh->nZ; ++iZ)
    {
      oneGroupXS->sigT(iZ,iR) = sigT(iZ,iR,0);      
      oneGroupXS->rSigTR(iZ,iR) = sigT(iZ,iR,0);      
      oneGroupXS->zSigTR(iZ,iR) = sigT(iZ,iR,0);      
      oneGroupXS->sigS(iZ,iR) = sigS(iZ,iR,0,0);      
      oneGroupXS->sigF(iZ,iR) = sigF(iZ,iR,0);      
      oneGroupXS->neutV(iZ,iR) = neutV(0);      
      oneGroupXS->rNeutV(iZ,iR) = neutV(0);      
      oneGroupXS->zNeutV(iZ,iR) = neutV(0);      
    }
  }

  for (int iR = 0; iR < mesh->nR; ++iR)
  {
    oneGroupXS->zSigTR(mesh->nZ,iR) = sigT(mesh->nZ-1,iR,0);      
    oneGroupXS->zNeutV(mesh->nZ,iR) = neutV(0);      
  }

  for (int iZ = 0; iZ < mesh->nZ; ++iZ)
  {
    oneGroupXS->rSigTR(iZ,mesh->nR) = sigT(iZ,mesh->nR-1,0);      
    oneGroupXS->rNeutV(iZ,mesh->nR) = neutV(0);      
  }

};
//==============================================================================

//==============================================================================
/// Print out material map and all entries in material bank

void Materials::edit()
{
  cout << "Material map: " << endl;
  cout << matMap << endl;
  cout << endl;
  for (int iCount = 0; iCount < matBank.size(); ++iCount){
    matBank[iCount]->edit();
  }

};
//==============================================================================
