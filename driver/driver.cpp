// test.cpp

#include <iostream>
#include "../libs/materials.h"
#include "../libs/multigroupQD.h"
#include "../libs/singlegroupQD.h"
#include "../libs/transport.h"
#include "../libs/mesh.h"

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
