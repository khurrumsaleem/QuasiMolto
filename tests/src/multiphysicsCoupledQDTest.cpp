
#include "../../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../../libs/Mesh.h"
#include "../../libs/Materials.h"
#include "../../libs/MultiGroupTransport.h"
#include "../../libs/MultiGroupQD.h"
#include "../../libs/MultiPhysicsCoupledQD.h"

using namespace std;

int main(int argc, char** argv)
{
  PetscInitialize(&argc,&argv,(char*)0,"Test");
  
  YAML::Node * input;
  input = new YAML::Node;
  *input = YAML::LoadFile("inputs/1-group-test.yaml");
  PetscErrorCode ierr;

  // initialize mesh object
  Mesh * myMesh;
  myMesh = new Mesh(input);

  // initialize materials object
  Materials * myMaterials;
  myMaterials = new Materials(myMesh,input);

  // initialize multigroup transport object
  MultiGroupTransport * myMGT; 
  myMGT = new MultiGroupTransport(myMaterials,myMesh,input);

  // initialize multigroup quasidiffusion object
  MultiGroupQD * myMGQD; 
  myMGQD = new MultiGroupQD(myMaterials,myMesh,input);

  // Initialize MultiPhysicsCoupledQD
  MultiPhysicsCoupledQD * myMPQD; 
  myMPQD = new MultiPhysicsCoupledQD(myMaterials,myMesh,input);

  cout << "Initialized multiphysics-coupled quasidiffusion class" << endl;
  
  ierr = PetscFinalize();
}
