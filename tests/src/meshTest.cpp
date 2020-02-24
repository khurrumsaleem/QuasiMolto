#include "../../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "../../libs/Mesh.h"

using namespace std;

int main()
{
  YAML::Node * input;
  input = new YAML::Node;
  *input = YAML::LoadFile("inputs/1-group-test.yaml");

  // initialize mesh object
  Mesh * myMesh;
  myMesh = new Mesh(input);
}
