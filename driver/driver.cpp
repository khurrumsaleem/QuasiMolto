// test.cpp

#include <iostream>
#include "../libs/Materials.h"
#include "../libs/MultiGroupQD.h"
#include "../libs/SingleGroupQD.h"
#include "../libs/Transport.h"
#include "../libs/Mesh.h"

using namespace std;

int main(void) {

     cout << "Print from driver" << endl;

     printMat();
     printMultiGroupQD();
     printSingleGroupQD();
     printTransport();

     Mesh myMesh; 
     myMesh.calcQuadSet();
     myMesh.printQuadSet();

     return(0);
}
