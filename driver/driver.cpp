// test.cpp

static char help[] = "Solves CFR reactor kinetics problems.\n\n";

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

enum SolveId 
{
  transient,
  steadyState,
  steadyStateThenTransient,
  steadyStateThenPseudoTransient,
  numSolveIds
};

const char* solveTypes[] =
{ 
  "transient",
  "steady_state",
  "steady_state_then_transient",
  "steady_state_then_pseudo_transient"
};

using namespace std;

// Declare various solver types
void transientSolve(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void steadyStateSolve(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void steadyStateThenTransientSolve(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void steadyStateThenPseudoTransientSolve(Materials * myMaterials,\
  Mesh * myMesh,\
  YAML::Node * input);
void printLogo();

int main(int argc, char** argv) {

  PetscMPIInt size;
  PetscErrorCode ierr; 
  string solveString;
  
  // Print logo
  printLogo();

  // Initialize PETSc
  PetscInitialize(&argc,&argv,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  // Get input file
  YAML::Node * input;
  input = new YAML::Node;
  if (argc>1)
    *input = YAML::LoadFile(argv[1]);
  else 
    *input = YAML::LoadFile("input.yaml");

  // Initialize mesh object
  Mesh * myMesh; 
  myMesh = new Mesh(input);
  if (myMesh->verbose)
    myMesh->printQuadSet();

  // Initialize materials object
  Materials * myMaterials;
  myMaterials = new Materials(myMesh,input);

  // Initialize multigroup transport object
  MultiGroupTransport * myMGT; 
  myMGT = new MultiGroupTransport(myMaterials,myMesh,input);

  // Initialize multigroup quasidiffusion object
  MultiGroupQD * myMGQD; 
  myMGQD = new MultiGroupQD(myMaterials,myMesh,input);

  // Initialize T2QD coupling object
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

  // Launch solver based on user-specified solve type
  if ((*input)["parameters"]["solve type"])
  {
    // Read in user-specified solve type  
    solveString=(*input)["parameters"]["solve type"].as<string>();
    cout << "User specified solve type: " << solveString << endl;

    // Transient solve
    if (solveString == solveTypes[transient])
      transientSolve(myMaterials,myMesh,input);
    // Steady state solve
    else if (solveString == solveTypes[steadyState])
      steadyStateSolve(myMaterials,myMesh,input);
    // Steady state initial solve followed by a transient solve
    else if (solveString == solveTypes[steadyStateThenTransient])
      steadyStateThenTransientSolve(myMaterials,myMesh,input);
    // Steady state initial solve followed by a pseudo transient solve 
    else if (solveString == solveTypes[steadyStateThenPseudoTransient])
      steadyStateThenPseudoTransientSolve(myMaterials,myMesh,input);
    // User-specified solve type is not supported
    else
    {
      cout << "Solve type " << solveString << " not recognized." << endl;
      cout << "Supported solve types are: " << endl; 
      for (int supportedSolveId = 0; 
           supportedSolveId != numSolveIds; 
           supportedSolveId++)
      {
        cout << "  ";
        cout << solveTypes[supportedSolveId] << endl;
      }
    }
  }
  // Default to steady state solve if solve type is not specified
  else
  {
    cout << "Defaulting to steady state solve." << endl;
    steadyStateSolve(myMaterials,myMesh,input);
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

//==============================================================================
/// Execute a transient solve
///
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] input YAML input object for the simulation
void transientSolve(Materials * myMaterials,\
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

  myMLCoupling->solveTransient();

  cout << "Completed multilevel transient solve." << endl;

}
//==============================================================================

//==============================================================================
/// Execute a steady state solve
///
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] input YAML input object for the simulation
void steadyStateSolve(Materials * myMaterials,\
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

  myMLCoupling->solveSteadyStateResidualBalance(true);

  cout << "Completed multilevel steady state solve." << endl;

  // Delete pointers
  delete myMGT;
  delete myMGQD;
  delete myMPQD;
  delete myMLCoupling;
}
//==============================================================================

//==============================================================================
/// Execute a steady state solve, and then use the steady state solution as the
/// the initial condition for a transient solve
///
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] input YAML input object for the simulation
void steadyStateThenTransientSolve(Materials * myMaterials,\
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

  cout << "Completed multilevel solve." << endl;

  // Delete pointers
  delete myMGT;
  delete myMGQD;
  delete myMPQD;
  delete myMLCoupling;
}
//==============================================================================

//==============================================================================
/// Execute a steady state solve, and then use the steady state solution as the
/// the initial condition for a pseudo transient solve. In the pseudo transient
/// solve, the precursor balance and heat transfer equations have time 
/// derivatives, but the neutron transport equations utilizes an eigenvalue
///
/// @param [in] myMaterials Materials object for the simulation
/// @param [in] myMesh Mesh object for the simulation
/// @param [in] input YAML input object for the simulation
void steadyStateThenPseudoTransientSolve(Materials * myMaterials,\
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

  myMLCoupling->solveSteadyStatePseudoTransient(true);

  cout << "Completed multilevel solve." << endl;

  // Delete pointers
  delete myMGT;
  delete myMGQD;
  delete myMPQD;
  delete myMLCoupling;
}
//==============================================================================

//==============================================================================
void printLogo()
{

  printf("  ____                          _   __  __           _   _           \n");
  printf(" / __ \\                        (_) |  \\/  |         | | | |          \n");
  printf("| |  | |  _   _    __ _   ___   _  | \\  / |   ___   | | | |_    ___  \n");
  printf("| |  | | | | | |  / _` | / __| | | | |\\/| |  / _ \\  | | | __|  / _ \\ \n");
  printf("| |__| | | |_| | | (_| | \\__ \\ | | | |  | | | (_) | | | | |_  | (_) |\n");
  printf(" \\___\\_\\  \\__,_|  \\__,_| |___/ |_| |_|  |_|  \\___/  |_|  \\__|  \\___/ \n");
  printf("                                                                     \n");
  printf("An open source research code for circulating fuel reactor kinetics \n\n");
  printf("BSD 3-Clause License\n");
  printf("Copyright (c) 2019, Aaron James Reynolds \n");
  printf("All rights reserved. \n\n");
}
//==============================================================================
