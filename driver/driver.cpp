// test.cpp

#include <iostream>
#include "../libs/Materials.h"
#include "../libs/MultiGroupQD.h"
#include "../libs/SingleGroupQD.h"
#include "../libs/Transport.h"
#include "../libs/Mesh.h"
#include "../libs/StartingAngle.h"
#include "../libs/Material.h"
#include "../libs/MultiGroupTransport.h"
#include "../libs/SingleGroupTransport.h"
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

     Materials * myMaterials;
     myMaterials = new Materials(myMesh,input);
    
     StartingAngle myStartingAngle(myMesh,myMaterials,input);
     myStartingAngle.calcStartingAngle();
     myMaterials->edit();

     MultiGroupTransport * myMGT; 
     myMGT = new MultiGroupTransport(myMaterials,myMesh,input);

     return(0);
}
