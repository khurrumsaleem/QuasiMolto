#include "../../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../../libs/Mesh.h"
#include "../../libs/Materials.h"
#include "../../libs/MultiGroupTransport.h"

using namespace std;

int main()
{
  YAML::Node * input;
  input = new YAML::Node;
  *input = YAML::LoadFile("inputs/1-group-test.yaml");

  // initialize mesh object
  Mesh * myMesh;
  myMesh = new Mesh(input);

  // initialize materials object
  Materials * myMaterials;
  myMaterials = new Materials(myMesh,input);

  // initialize multigroup transport object
  MultiGroupTransport * myMGT;
  myMGT = new MultiGroupTransport(myMaterials,myMesh,input);
}
