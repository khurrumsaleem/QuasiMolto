// test.cpp

static char help[] = "Solves CFR reactor kinetics problem with KSP.\n\n";

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
#include "../libs/MultiGroupQDToMultiPhysicsQDCoupling.h"
#include "../libs/MultilevelCoupling.h"
#include "../libs/PETScWrapper.h"
#include "../libs/MMS.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"

using namespace std;
void testHeatTransfer(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void testMultiGroupPrecursor(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void testMultiPhysicsCoupledQD(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void testMultiGroupToGreyGroupCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void testMultilevelCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void testSteadyState(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void testSteadyStateThenTransient(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void testPETScCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);

int main(int argc, char** argv) {

  PetscMPIInt size;
  PetscErrorCode ierr; 
  string solveType;

  // initialize PETSc
  PetscInitialize(&argc,&argv,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

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

  // Set number of procs
  if ((*input)["parameters"]["nprocs"]){
    int nProcs = (*input)["parameters"]["nprocs"].as<int>();
    Eigen::initParallel();
    omp_set_num_threads(nProcs);
    Eigen::setNbThreads(nProcs);
  }
  else
  {
    omp_set_num_threads(1);
    Eigen::setNbThreads(1);
  }

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
    else if (solveType == "testMultiPhysicsCoupledQD")
      testMultiPhysicsCoupledQD(myMaterials,myMesh,input);
    else if (solveType == "testMultiGroupToGreyGroupCoupling")
      testMultiGroupToGreyGroupCoupling(myMaterials,myMesh,input);
    else if (solveType == "testMultilevelCoupling")
      testMultilevelCoupling(myMaterials,myMesh,input);
    else if (solveType == "testSteadyState")
      testSteadyState(myMaterials,myMesh,input);
    else if (solveType == "testSteadyStateThenTransient")
      testSteadyStateThenTransient(myMaterials,myMesh,input);
    else if (solveType == "testPETScCoupling")
      testPETScCoupling(myMaterials,myMesh,input);
    else
      myMGT->solveTransportOnly();
  }
  else
  {
    myMGT->solveTransportOnly();
  }

  // Delete pointers

  delete myMesh;
  delete myMaterials;
  delete myMGT; 
  delete myMGQD; 
  delete myT2QD; 
  delete myMMS;

  // Finalize PETSc
  ierr = PetscFinalize();

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

  myMPQD->mgdnp->buildRecircLinearSystem();
  cout << "recirculation A: " << endl;
  cout << myMPQD->mgdnp->recircA << endl;
  myMPQD->mgdnp->solveRecircLinearSystem();

  cout << myMPQD->mgdnp->beta << endl;
  cout << "recirculation Z: " << myMesh->recircZ << endl;
  
  myMPQD->mgdnp->DNPs[0]->assignBoundaryIndices();
  myMPQD->mgdnp->DNPs[0]->updateBoundaryConditions();
  cout << "Updated boundary conditions" << endl;
 
  myMPQD->mgdnp->DNPs[0]->calcRecircDNPFluxes();
  myMPQD->mgdnp->DNPs[0]->calcCoreDNPFluxes();
  cout << "Setting size of A and b..." << endl;
  myMPQD->A.resize(myMesh->nZ*myMesh->nR,myMesh->nZ*myMesh->nR);
  myMPQD->b.resize(myMesh->nZ*myMesh->nR);
  cout << "Set size of A and b" << endl;
  myMPQD->mgdnp->DNPs[0]->buildCoreLinearSystem();
  cout << "built core system" << endl;
  cout << "A" << endl;
  cout << myMPQD->A << endl;
  cout << "b" << endl;
  cout << myMPQD->b << endl;
  myMPQD->solveLinearSystem();
  cout << "x: " << endl;   
  cout << myMPQD->x << endl;
  cout << "building recirc system..." << endl;
  myMPQD->mgdnp->recircA.resize(myMesh->nZrecirc*myMesh->nR,\
    myMesh->nZrecirc*myMesh->nR);
  myMPQD->mgdnp->recircb.resize(myMesh->nZrecirc*myMesh->nR);
  myMPQD->mgdnp->DNPs[0]->buildRecircLinearSystem();
  cout << "built recirc system" << endl;
  cout << "recircA:" << endl;
  cout << myMPQD->mgdnp->recircA << endl;
  cout << "recircb:" << endl;
  cout << myMPQD->mgdnp->recircb << endl;
  myMPQD->mgdnp->solveRecircLinearSystem();
  cout << "recircx: " << endl;   
  cout << myMPQD->mgdnp->recircx << endl;
}

void testMultiPhysicsCoupledQD(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input){
  
  MultiPhysicsCoupledQD * myMPQD; 
  myMPQD = new MultiPhysicsCoupledQD(myMaterials,myMesh,input);
  myMPQD->solveTransient();
}

void testMultiGroupToGreyGroupCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input){
  
  // initialize multigroup quasidiffusion object
  MultiGroupQD * myMGQD; 
  myMGQD = new MultiGroupQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MultiPhysicsCoupledQD * myMPQD; 
  myMPQD = new MultiPhysicsCoupledQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MGQDToMPQDCoupling * myMGToGG; 
  myMGToGG = new MGQDToMPQDCoupling(myMesh,myMaterials,input,myMPQD,myMGQD);

  cout << "Initialized MGQDToMGQDCoupling" << endl;

  // Collapse nuclear data
  myMGToGG->solveTransient();
  myMaterials->oneGroupXS->print();
  myMPQD->ggqd->printBCParams();
  
  cout << "Collapsed nuclear data" << endl;
}

void testMultilevelCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input){
  
  // initialize multigroup transport object
  MultiGroupTransport * myMGT; 
  myMGT = new MultiGroupTransport(myMaterials,myMesh,input);

  // initialize multigroup quasidiffusion object
  MultiGroupQD * myMGQD; 
  myMGQD = new MultiGroupQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MultiPhysicsCoupledQD * myMPQD; 
  myMPQD = new MultiPhysicsCoupledQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MultilevelCoupling * myMLCoupling; 
  myMLCoupling = new MultilevelCoupling(myMesh,myMaterials,input,myMGT,myMGQD,\
      myMPQD);

  cout << "Initialized multilevel coupling" << endl;

  myMLCoupling->solveTransient();
  
  cout << "Completed multilevel solve" << endl;
  
//  myMaterials->oneGroupXS->print();
//  myMPQD->ggqd->printBCParams();
}

void testSteadyState(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input){
  
  // initialize multigroup transport object
  MultiGroupTransport * myMGT; 
  myMGT = new MultiGroupTransport(myMaterials,myMesh,input);

  // initialize multigroup quasidiffusion object
  MultiGroupQD * myMGQD; 
  myMGQD = new MultiGroupQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MultiPhysicsCoupledQD * myMPQD; 
  myMPQD = new MultiPhysicsCoupledQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MultilevelCoupling * myMLCoupling; 
  myMLCoupling = new MultilevelCoupling(myMesh,myMaterials,input,myMGT,myMGQD,\
      myMPQD);

  // Set state to zero for initial steady state solve
  myMesh->state=0;
  
  cout << "Initialized steady state solve" << endl;

  myMLCoupling->solveSteadyStateResidualBalance(true);

  cout << "Completed steady state solve" << endl;
  
  // Delete pointers

  delete myMGT;
  delete myMGQD;
  delete myMPQD;
  delete myMLCoupling;

}

void testSteadyStateThenTransient(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input){
  
  // initialize multigroup transport object
  MultiGroupTransport * myMGT; 
  myMGT = new MultiGroupTransport(myMaterials,myMesh,input);

  // initialize multigroup quasidiffusion object
  MultiGroupQD * myMGQD; 
  myMGQD = new MultiGroupQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MultiPhysicsCoupledQD * myMPQD; 
  myMPQD = new MultiPhysicsCoupledQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MultilevelCoupling * myMLCoupling; 
  myMLCoupling = new MultilevelCoupling(myMesh,myMaterials,input,myMGT,myMGQD,\
      myMPQD);

  cout << "Starting solves..." << endl;

  myMLCoupling->solveSteadyStateTransientResidualBalance(true);
  
  cout << "Completed multilevel solve" << endl;

  // Delete pointers

  delete myMGT;
  delete myMGQD;
  delete myMPQD;
  delete myMLCoupling;

}

void testPETScCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input){
  
  // initialize multigroup transport object
  MultiGroupTransport * myMGT; 
  myMGT = new MultiGroupTransport(myMaterials,myMesh,input);

  // initialize multigroup quasidiffusion object
  MultiGroupQD * myMGQD; 
  myMGQD = new MultiGroupQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MultiPhysicsCoupledQD * myMPQD; 
  myMPQD = new MultiPhysicsCoupledQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MultilevelCoupling * myMLCoupling; 
  myMLCoupling = new MultilevelCoupling(myMesh,myMaterials,input,myMGT,myMGQD,\
      myMPQD);

  cout << "Build PETSc MGQD linear system...";
  myMGQD->buildSteadyStateLinearSystem_p();
  cout << " done." << endl;
  
  cout << "Solve system...";
  myMGQD->solveLinearSystem_p();
  cout << " done." << endl;
  
  cout << "Build Eigen MGQD linear system...";
  myMGQD->buildSteadyStateLinearSystem();
  cout << " done." << endl;
  
  cout << "Solve system...";
  myMGQD->QDSolve->solveSuperLU();
  cout << " done." << endl;


  // Delete pointers

  delete myMGT;
  delete myMGQD;
  delete myMPQD;
  delete myMLCoupling;

}
