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
int testMGQDPETScCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
int testELOTPETScCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
int testSteadyStateMultilevelPETScCoupling(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void testTransientMultilevelPETScCoupling(Materials * myMaterials,\
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
    else if (solveType == "transient")
      testMultilevelCoupling(myMaterials,myMesh,input);
    else if (solveType == "steady_state")
      testSteadyState(myMaterials,myMesh,input);
    else if (solveType == "testSteadyStateThenTransient")
      testSteadyStateThenTransient(myMaterials,myMesh,input);
    else if (solveType == "testMGQDPETScCoupling")
      testMGQDPETScCoupling(myMaterials,myMesh,input);
    else if (solveType == "testELOTPETScCoupling")
      testELOTPETScCoupling(myMaterials,myMesh,input);
    else if (solveType == "testSteadyStateMultilevelPETScCoupling")
      testSteadyStateMultilevelPETScCoupling(myMaterials,myMesh,input);
    else if (solveType == "testTransientMultilevelPETScCoupling")
      testTransientMultilevelPETScCoupling(myMaterials,myMesh,input);
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

  cout << "Initialized multilevel transient solve." << endl;

  if (myMesh->petsc)
    myMLCoupling->solveTransient_p();
  else
    myMLCoupling->solveTransient();

  cout << "Completed multilevel transient solve." << endl;

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

  cout << "Initialized multilevel steady state solve." << endl;

  if (myMesh->petsc)
    myMLCoupling->solveSteadyStateResidualBalance_p(true);
  else
    myMLCoupling->solveSteadyStateResidualBalance(true);


  cout << "Completed multilevel steady state solve." << endl;

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

int testMGQDPETScCoupling(Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * input){

  PetscErrorCode ierr;

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
  //
  //  /* PETSc steady state*/  
  //  cout << "Build PETSc MGQD linear system...";
  //  myMGQD->buildSteadyStateLinearSystem_p();
  //  cout << " done." << endl;
  //  
  //  cout << "Solve system...";
  //  myMGQD->solveLinearSystem_p();
  //  cout << " done." << endl;
  //
  //  cout << "Build system to back calculate current...";
  //  myMGQD->buildSteadyStateBackCalcSystem_p();
  //  cout << " done." << endl;
  //  
  //  cout << "Solve system...";
  //  myMGQD->backCalculateCurrent_p();
  //  cout << " done." << endl;
  //  myMGQD->getFluxes();
  //  myMGQD->writeVars();
  //
  //  /* EIGEN steady state*/  
  //  cout << "Build Eigen MGQD linear system...";
  //  myMGQD->buildSteadyStateLinearSystem();
  //  cout << " done." << endl;
  //  
  //  cout << "Solve system...";
  //  myMGQD->QDSolve->solveSuperLU();
  //  cout << " done." << endl;
  //  
  //  cout << "Build system to back calculate current...";
  //  myMGQD->buildSteadyStateBackCalcSystem();
  //  cout << " done." << endl;
  //  
  //  cout << "Solve system...";
  //  myMGQD->backCalculateCurrent();
  //  cout << " done." << endl;
  //
  //  /* PETSc transient*/  
  //  cout << "Build transient PETSc MGQD linear system...";
  //  myMGQD->buildLinearSystem_p();
  //  cout << " done." << endl;
  //  
  //  cout << "Solve system...";
  //  myMGQD->solveLinearSystem_p();
  //  cout << " done." << endl;
  //  //ierr = VecView(myMGQD->QDSolve->x_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //
  //  cout << "Build system to back calculate current...";
  //  myMGQD->buildBackCalcSystem_p();
  //  cout << " done." << endl;
  //
  //  cout << "Solve system...";
  //  myMGQD->backCalculateCurrent_p();
  //  cout << " done." << endl;
  //  ierr = VecView(myMGQD->QDSolve->currPast_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //  //myMGQD->getFluxes();
  //  //myMGQD->writeVars();
  //
  //  /* EIGEN transient*/  
  //  cout << "Build Eigen MGQD linear system...";
  //  myMGQD->buildLinearSystem();
  //  cout << " done." << endl;
  //  
  //  cout << "Solve system...";
  //  myMGQD->QDSolve->solveSuperLU();
  //  cout << " done." << endl;
  //  //cout << myMGQD->QDSolve->x << endl;
  //  
  //  cout << "Build system to back calculate current...";
  //  myMGQD->buildBackCalcSystem();
  //  cout << " done." << endl;
  //  
  //  cout << "Solve system...";
  //  myMGQD->backCalculateCurrent();
  //  cout << " done." << endl;
  //  cout << myMGQD->QDSolve->currPast<< endl;

  /* PETSc transient*/  
  if (myMesh->petsc)
  {
    cout << "Run PETSc transient...";
    myMGQD->solveMGQDOnly_p();
    cout << " done." << endl;
  }
  else
  {
    cout << "Run Eigen transient...";
    myMGQD->solveMGQDOnly();
    cout << " done." << endl;
  }

  // Delete pointers
  delete myMGT;
  delete myMGQD;
  delete myMPQD;
  delete myMLCoupling;

}

int testELOTPETScCoupling(Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * input){

  PetscErrorCode ierr;

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

  // Test GGSolver setFlux
  //ierr = VecView(myMPQD->xPast_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /* STEADY STATE */

  // Test GGSolver steadyStateCurrent functions
  double testDub = 0.1;
  int testInt = 0;
  //myMPQD->ggqd->GGSolver->steadyStateSouthCurrent_p(testDub,\
  //    testInt,testInt,testInt);
  //cout << "south current" << endl;
  //myMPQD->ggqd->GGSolver->steadyStateNorthCurrent_p(testDub,\
  //    testInt,testInt,testInt);
  //cout << "north current" << endl;
  //myMPQD->ggqd->GGSolver->steadyStateWestCurrent_p(testDub,\
  //    testInt,testInt,testInt);
  //cout << "west current" << endl;
  //myMPQD->ggqd->GGSolver->steadyStateEastCurrent_p(testDub,\
  //    testInt,testInt,testInt);
  //cout << "east current" << endl;

  //// Test zeroth moment, radial, and axial current conditions
  //myMPQD->ggqd->GGSolver->assertSteadyStateZerothMoment_p(testInt,testInt,testInt);
  //cout << "zeroth moment" << endl;
  //myMPQD->ggqd->GGSolver->applySteadyStateRadialBoundary_p(testInt,testInt,testInt);
  //cout << "radial boundary" << endl;
  //myMPQD->ggqd->GGSolver->applySteadyStateAxialBoundary_p(testInt,testInt,testInt);
  //cout << "axial boundary" << endl;

  //// Test steady-state boundary conditions on current
  //myMPQD->ggqd->GGSolver->assertSteadyStateWCurrentBC_p(testInt,testInt,testInt);
  //cout << "west current bc" << endl;
  //myMPQD->ggqd->GGSolver->assertSteadyStateECurrentBC_p(testInt,testInt,testInt);
  //cout << "east current bc" << endl;
  //myMPQD->ggqd->GGSolver->assertSteadyStateNCurrentBC_p(testInt,testInt,testInt);
  //cout << "north current bc" << endl;
  //myMPQD->ggqd->GGSolver->assertSteadyStateSCurrentBC_p(testInt,testInt,testInt);
  //cout << "south current bc" << endl;

  //// Test steady-state boundary conditions on current
  //myMPQD->ggqd->GGSolver->assertSteadyStateNGoldinBC_p(testInt,testInt,testInt);
  //cout << "north current goldin bc" << endl;
  //myMPQD->ggqd->GGSolver->assertSteadyStateSGoldinBC_p(testInt,testInt,testInt);
  //cout << "south current goldin bc" << endl;
  //myMPQD->ggqd->GGSolver->assertSteadyStateEGoldinBC_p(testInt,testInt,testInt);
  //cout << "east current goldin bc" << endl;

  //// Test steady-state P1 boundary conditions on current
  //myMPQD->ggqd->GGSolver->assertSteadyStateNGoldinP1BC_p(testInt,testInt,testInt);
  //cout << "north current goldin p1 bc" << endl;
  //myMPQD->ggqd->GGSolver->assertSteadyStateSGoldinP1BC_p(testInt,testInt,testInt);
  //cout << "south current goldin p1 bc" << endl;
  //myMPQD->ggqd->GGSolver->assertSteadyStateEGoldinP1BC_p(testInt,testInt,testInt);
  //cout << "east current goldin p1 bc" << endl;

  //// Test steady-state P1 boundary conditions on current
  //myMPQD->ggqd->GGSolver->assertNFluxBC_p(testInt,testInt,testInt);
  //cout << "north flux bc" << endl;
  //myMPQD->ggqd->GGSolver->assertSFluxBC_p(testInt,testInt,testInt);
  //cout << "south flux bc" << endl;
  //myMPQD->ggqd->GGSolver->assertEFluxBC_p(testInt,testInt,testInt);
  //cout << "east flux bc" << endl;
  //myMPQD->ggqd->GGSolver->assertWFluxBC_p(testInt,testInt,testInt);
  //cout << "west flux bc" << endl;

  //// Test build steady state QD linear system call
  //myMPQD->ggqd->GGSolver->formSteadyStateLinearSystem_p();
  //cout << "form ggqd steady state linear system" << endl;

  //// Test build steady state heat transfer linear system call
  //myMPQD->heat->buildSteadyStateLinearSystem_p();
  //cout << "form heat transfer steady state linear system" << endl;

  //// Test build steady state precursor balance linear system call
  //myMPQD->mgdnp->buildSteadyStateCoreLinearSystem_p();
  //cout << "form mgdnp transfer steady state linear system" << endl;

  // Test calcSteadyStateCurrent functions
  //myMPQD->ggqd->GGSolver->calcSteadyStateSouthCurrent_p(testInt,testInt,testInt);
  //cout << "calc south current" << endl;
  //myMPQD->ggqd->GGSolver->calcSteadyStateNorthCurrent_p(testInt,testInt,testInt);
  //cout << "calc north current" << endl;
  //myMPQD->ggqd->GGSolver->calcSteadyStateWestCurrent_p(testInt,testInt,testInt);
  //cout << "calc west current" << endl;
  //myMPQD->ggqd->GGSolver->calcSteadyStateEastCurrent_p(testInt,testInt,testInt);
  //cout << "calc east current" << endl;

  // Build full ELOT system
  //myMPQD->buildSteadyStateLinearSystem_p();
  ////ierr = VecView(myMPQD->b_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //myMPQD->solve_p();
  //ierr = VecView(myMPQD->x_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //myMPQD->updateSteadyStateVarsAfterConvergence_p();
  if (myMesh->petsc)
  {
    cout << "Run ELOT PETSc steady state solve...";
    myMPQD->solveSteadyState_p();
    //ierr = VecView(myMPQD->x_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    //ierr = VecView(myMPQD->b_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    //ierr = MatView(myMPQD->A_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    //ierr = MatView(myMPQD->ggqd->GGSolver->C_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    //cout << " done." << endl;
  }
  else
  {
    cout << "Run ELOT eigen steady state solve...";
    myMPQD->solveSteadyState();
    //cout << "x" << endl;
    //cout << myMPQD->x << endl;
    //cout << "b" << endl;
    //cout << myMPQD->b << endl;
    //cout << "A" << endl;
    //cout << myMPQD->A << endl;
    //cout << "C" << endl;
    //cout << myMPQD->ggqd->GGSolver->C << endl;
    //cout << " done." << endl;
  }

  /* TRANSIENT */

  //VecScatter     ctx;
  //VecScatterCreateToAll(myMPQD->ggqd->GGSolver->currPast_p,&ctx,&(myMPQD->ggqd->GGSolver->currPast_p_seq));
  //VecScatterBegin(ctx,myMPQD->ggqd->GGSolver->currPast_p,myMPQD->ggqd->GGSolver->currPast_p_seq,\
  //    INSERT_VALUES,SCATTER_FORWARD);
  //VecScatterEnd(ctx,myMPQD->ggqd->GGSolver->currPast_p,myMPQD->ggqd->GGSolver->currPast_p_seq,\
  INSERT_VALUES,SCATTER_FORWARD);

  //myMPQD->ggqd->GGSolver->southCurrent_p(testDub,\
  //    testInt,testInt,testInt);
  //cout << "south current" << endl;
  //myMPQD->ggqd->GGSolver->northCurrent_p(testDub,\
  //    testInt,testInt,testInt);
  //cout << "north current" << endl;
  //myMPQD->ggqd->GGSolver->westCurrent_p(testDub,\
  //    testInt,testInt,testInt);
  //cout << "west current" << endl;
  //myMPQD->ggqd->GGSolver->eastCurrent_p(testDub,\
  //    testInt,testInt,testInt);
  //cout << "east current" << endl;

  //// Test steady-state boundary conditions on current
  //myMPQD->ggqd->GGSolver->assertWCurrentBC_p(testInt,testInt,testInt);
  //cout << "west current bc" << endl;
  //myMPQD->ggqd->GGSolver->assertECurrentBC_p(testInt,testInt,testInt);
  //cout << "east current bc" << endl;
  //myMPQD->ggqd->GGSolver->assertNCurrentBC_p(testInt,testInt,testInt);
  //cout << "north current bc" << endl;
  //myMPQD->ggqd->GGSolver->assertSCurrentBC_p(testInt,testInt,testInt);
  //cout << "south current bc" << endl;

  //// Test steady-state boundary conditions on current
  //myMPQD->ggqd->GGSolver->assertNGoldinBC_p(testInt,testInt,testInt);
  //cout << "north current goldin bc" << endl;
  //myMPQD->ggqd->GGSolver->assertSGoldinBC_p(testInt,testInt,testInt);
  //cout << "south current goldin bc" << endl;
  //myMPQD->ggqd->GGSolver->assertEGoldinBC_p(testInt,testInt,testInt);
  //cout << "east current goldin bc" << endl;

  //// Test steady-state P1 boundary conditions on current
  //myMPQD->ggqd->GGSolver->assertNGoldinP1BC_p(testInt,testInt,testInt);
  //cout << "north current goldin p1 bc" << endl;
  //myMPQD->ggqd->GGSolver->assertSGoldinP1BC_p(testInt,testInt,testInt);
  //cout << "south current goldin p1 bc" << endl;
  //myMPQD->ggqd->GGSolver->assertEGoldinP1BC_p(testInt,testInt,testInt);
  //cout << "east current goldin p1 bc" << endl;

  //myMPQD->ggqd->GGSolver->assertZerothMoment_p(testInt,testInt,testInt);
  //cout << "zeroth moment" << endl;
  //myMPQD->ggqd->GGSolver->applyRadialBoundary_p(testInt,testInt,testInt);
  //cout << "radial boundary" << endl;
  //myMPQD->ggqd->GGSolver->applyAxialBoundary_p(testInt,testInt,testInt);
  //cout << "axial boundary" << endl;
  //
  //// Test build steady state QD linear system call
  //myMPQD->ggqd->GGSolver->formLinearSystem_p();
  //cout << "form ggqd linear system" << endl;
  //
  //// Test build heat transfer linear system call
  //myMPQD->heat->buildLinearSystem_p();
  //cout << "form heat transfer linear system" << endl;

  //// Test build steady state precursor balance linear system call
  //myMPQD->mgdnp->buildCoreLinearSystem_p();
  //cout << "form mgdnp core linear system" << endl;

  // Test build steady state precursor balance linear system call
  //myMPQD->mgdnp->buildRecircLinearSystem_p();
  //cout << "form mgdnp recirc linear system" << endl;
  //myMPQD->mgdnp->solveRecircLinearSystem_p();
  //cout << "solved mgdnp recirc linear system" << endl;

  // Test calcSteadyStateCurrent functions
  //myMPQD->ggqd->GGSolver->calcSteadyStateSouthCurrent_p(testInt,testInt,testInt);
  //cout << "calc south current" << endl;
  //myMPQD->ggqd->GGSolver->calcSteadyStateNorthCurrent_p(testInt,testInt,testInt);
  //cout << "calc north current" << endl;
  //myMPQD->ggqd->GGSolver->calcSteadyStateWestCurrent_p(testInt,testInt,testInt);
  //cout << "calc west current" << endl;
  //myMPQD->ggqd->GGSolver->calcEastCurrent_p(testInt,testInt,testInt);
  //cout << "calc east current" << endl;

  //// Test formBackCalc functions
  //myMPQD->ggqd->GGSolver->formBackCalcSystem_p();
  //cout << "form back calc system" << endl;

  //PETSc system 
  //if (myMesh->petsc)
  //{
  //  myMPQD->solveTransient_p();
  //  //cout << "x_p" << endl;
  //  //ierr = VecView(myMPQD->x_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //  //cout << "b_p" << endl;
  //  //ierr = VecView(myMPQD->b_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //  //cout << "A_p" << endl;
  //  //ierr = MatView(myMPQD->A_p,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //}
  //// Eigen system 
  //else
  //{
  //  myMPQD->solveTransient();
  //  //cout << "x" << endl;
  //  //cout << myMPQD->x << endl;
  //  //cout << "b" << endl;
  //  //cout << myMPQD->b << endl;
  //  //cout << "A" << endl;
  //  //cout << myMPQD->A << endl;
  //  //cout << "C" << endl;
  //  //cout << myMPQD->ggqd->GGSolver->C << endl;
  //  //cout << " done." << endl;

  //}


  // Delete pointers
  delete myMGT;
  delete myMGQD;
  delete myMPQD;
  delete myMLCoupling;

}

int testSteadyStateMultilevelPETScCoupling(Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * input){

  PetscErrorCode ierr;

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

  if (myMesh->petsc)
    myMLCoupling->solveSteadyStateResidualBalance_p(true);
  else
    myMLCoupling->solveSteadyStateResidualBalance(true);

  cout << "Completed steady state solve" << endl;

  // Delete pointers
  delete myMGT;
  delete myMGQD;
  delete myMPQD;
  delete myMLCoupling;

}

void testTransientMultilevelPETScCoupling(Materials * myMaterials,\
    Mesh * myMesh,\
    YAML::Node * input){

  PetscErrorCode ierr;

  // initialize multigroup transport object
  MultiGroupTransport * myMGT; 
  myMGT = new MultiGroupTransport(myMaterials,myMesh,input);

  // initialize multigroup quasidiffusion object
  MultiGroupQD * myMGQD; 
  myMGQD = new MultiGroupQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MultiPhysicsCoupledQD * myMPQD; 
  myMPQD = new MultiPhysicsCoupledQD(myMaterials,myMesh,input);

  // Initialize MultilevelCoupling
  MultilevelCoupling * myMLCoupling; 
  myMLCoupling = new MultilevelCoupling(myMesh,myMaterials,input,myMGT,myMGQD,\
      myMPQD);

  if (myMesh->petsc)
    myMLCoupling->solveTransient_p();
  else
    myMLCoupling->solveTransient();

  cout << "Completed transient solve" << endl;

  // Delete pointers
  delete myMGT;
  delete myMGQD;
  delete myMPQD;
  delete myMLCoupling;

}
