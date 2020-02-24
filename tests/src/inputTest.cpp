#include "../../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"

using namespace std;

int main()
{
  YAML::Node * input;
  input = new YAML::Node;
  *input = YAML::LoadFile("inputs/1-group-test.yaml");
}
