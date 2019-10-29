#include "../TPLs/yaml-cpp/include/yaml-cpp/yaml.h"
#include "Mesh.h"

using namespace std; 
using namespace arma;

//==============================================================================
//! StartingAngle class that solves RZ neutron transport at the starting angles

class StartingAngle
{
        public:
        // public functions
        StartingAngle(Mesh * myMesh,YAML::Node * myInput);
	void calcStartingAngle();
        mat calckR(double myGamma);
        mat calckZ(double myGamma);
        mat calclR(double myGamma);
        mat calclZ(double myGamma);
        mat calct1(double myGamma);
        mat calct2(double myGamma);
        rowvec calcSubCellVol(int myiZ, int myiR);

        private:
        // private functions
        YAML::Node * input;
        Mesh * mesh;
};

//==============================================================================
