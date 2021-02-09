#include "../../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../../libs/Mesh.h"
#include "../../libs/Materials.h"
#include "../../libs/QuasidiffusionSolver.h"

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

  // initialize quasidiffusionsolver
  QDSolver * myQDSolver;
  myQDSolver = new QDSolver(myMesh,myMaterials,input);
  
  ierr = PetscFinalize();
}
