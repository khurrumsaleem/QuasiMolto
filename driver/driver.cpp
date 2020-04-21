// test.cpp

#include <iostream>
#include "../libs/Materials.h"
#include "../libs/MultiGroupQD.h"
#include "../libs/SingleGroupQD.h"
#include "../libs/Mesh.h"
#include "../libs/StartingAngle.h"
#include "../libs/Material.h"
#include "../libs/MultiGroupTransport.h"
#include "../libs/SingleGroupTransport.h"
#include "../libs/SimpleCornerBalance.h"
#include "../libs/QuasidiffusionSolver.h"
#include "../libs/TransportToQDCoupling.h"
#include "../libs/HeatTransfer.h"
#include "../libs/MultiGroupDNP.h"
#include "../libs/SingleGroupDNP.h"
#include "../libs/MultiPhysicsCoupledQD.h"
#include "../libs/MMS.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"

using namespace std;
void testHeatTransfer(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void testMultiGroupPrecursor(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);

int main(int argc, char** argv) {

  string solveType;

  // get input file
  YAML::Node * input;
  input = new YAML::Node;
  if (argc>1){
    *input = YAML::LoadFile(argv[1]);
  } else {
    *input = YAML::LoadFile("input.yaml");
  }
       
  // initialize mesh object
  Mesh * myMesh; 
  myMesh = new Mesh(input);
  myMesh->printQuadSet();

  // initialize materials object
  Materials * myMaterials;
  myMaterials = new Materials(myMesh,input);
      
  // initialize multigroup transport object
  MultiGroupTransport * myMGT; 
  myMGT = new MultiGroupTransport(myMaterials,myMesh,input);
  
  // initialize multigroup quasidiffusion object
  MultiGroupQD * myMGQD; 
  myMGQD = new MultiGroupQD(myMaterials,myMesh,input);
  
  // initialize T2QD coupling object
  TransportToQDCoupling * myT2QD; 
  myT2QD = new TransportToQDCoupling(myMaterials,myMesh,input,myMGT,myMGQD);

  MMS * myMMS;
  myMMS = new MMS(myMGT,myMesh,myMaterials,input);

  if ((*input)["parameters"]["solve type"]){

    solveType=(*input)["parameters"]["solve type"].as<string>();

    cout << solveType << endl;
  
    if (solveType == "MMS" or solveType == "mms")    
      myMMS->timeDependent();
    else if (solveType == "MGQD" or solveType == "mgqd") 
      myMGQD->solveMGQDOnly();
    else if (solveType == "TQD" or solveType == "TQD")
    { 
     // myMGT->solveTransportOnly();
     // myT2QD->calcEddingtonFactors();
     // myT2QD->calcBCs();
     // myMGQD->solveMGQDOnly();
     myT2QD->solveTransportWithQDAcceleration();
    }
    else if (solveType == "testHeatTransfer")
      testHeatTransfer(myMaterials,myMesh,input);
    else if (solveType == "testMultiGroupPrecursor")
      testMultiGroupPrecursor(myMaterials,myMesh,input);
    else
      myMGT->solveTransportOnly();
  }
  else
  {
    myMGT->solveTransportOnly();
  }
return(0);
}

void testHeatTransfer(Materials * myMaterials,Mesh * myMesh,YAML::Node * input){
  
  MultiPhysicsCoupledQD * myMPQD; 
  myMPQD = new MultiPhysicsCoupledQD(myMaterials,myMesh,input);
  myMPQD->heat->updateBoundaryConditions();
  myMPQD->heat->calcDiracs();
  cout << "diracs: " << endl;   
  cout << myMPQD->heat->dirac << endl;
  cout << "flux: " << endl;   
  cout << myMPQD->heat->flux << endl;
  myMPQD->heat->calcFluxes();
  myMaterials->oneGroupXS->sigF.setOnes(myMesh->nZ,myMesh->nR);
  myMPQD->A.resize(myMesh->nZ*myMesh->nR,myMesh->nZ*myMesh->nR);
  cout << "Set size of A." << endl;
  myMPQD->b.resize(myMesh->nZ*myMesh->nR);
  myMPQD->x.resize(myMesh->nZ*myMesh->nR);
  cout << "Set size of b." << endl;
  myMPQD->heat->buildLinearSystem();
  cout << "A: " << endl;   
  cout << myMPQD->A << endl;;
  cout << "b: " << endl;   
  cout << myMPQD->b << endl;;
  myMPQD->solveLinearSystem();
  cout << "x: " << endl;   
  cout << myMPQD->x << endl;;
  myMPQD->heat->getTemp();
  cout << myMPQD->heat->temp << endl;;
  
}

void testMultiGroupPrecursor(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input){
  
  MultiPhysicsCoupledQD * myMPQD; 
  MultiGroupDNP * myMGP; 
  myMPQD = new MultiPhysicsCoupledQD(myMaterials,myMesh,input);
  cout << myMPQD->dnps->beta << endl;
  cout << "recirculation Z: " << myMesh->recircZ << endl;
  myMPQD->dnps->DNPs[0]->assignBoundaryIndices();
  myMPQD->dnps->DNPs[0]->updateBoundaryConditions();
  cout << "Updated boundary conditions" << endl;
}
