// test.cpp

#include <iostream>
#include "../libs/Materials.h"
#include "../libs/MultiGroupQD.h"
#include "../libs/SingleGroupQD.h"
#include "../libs/Transport.h"
#include "../libs/Mesh.h"
#include "../libs/StartingAngle.h"
#include "../libs/Material.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"

using namespace std;

int main(void) {

     YAML::Node * input;
     input = new YAML::Node;
     *input = YAML::LoadFile("input.yaml");
     cout << "Print from driver" << endl;

     printMultiGroupQD();
     printSingleGroupQD();
     printTransport();

     Mesh * myMesh; 
     myMesh = new Mesh(input);
     myMesh->printQuadSet();
     StartingAngle myStartingAngle(myMesh,input);
     myStartingAngle.calcStartingAngle();

     Materials * myMaterials;
     myMaterials = new Materials(myMesh,input);
    
     return(0);
}
