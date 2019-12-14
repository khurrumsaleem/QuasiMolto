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
#include "../libs/SimpleCornerBalance.h"
#include "../libs/MMS.h"
#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"

using namespace std;

int main(int argc, char** argv) {

  string solveType;

  // get input file
  YAML::Node * input;
  input = new YAML::Node;
  if (argc>1){
    *input = YAML::LoadFile(argv[1]);
  } else {
    *input = YAML::LoadFile("input.yaml");
  }
       
  printMultiGroupQD();
  printSingleGroupQD();
  printTransport();

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

  MMS * myMMS;
  myMMS = new MMS(myMGT,myMesh,myMaterials,input);

  if ((*input)["parameters"]["solve type"]){

    solveType=(*input)["parameters"]["solve type"].as<string>();

    cout << solveType << endl;
  
    if (solveType == "MMS" or solveType == "mms")    
      myMMS->timeDependent();
    else 
      myMGT->solveTransportOnly();

  }
  else
    myMGT->solveTransportOnly();

return(0);
}
