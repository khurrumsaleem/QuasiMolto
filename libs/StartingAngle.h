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

        private:
        // private functions
        YAML::Node * input;
        Mesh * mesh;
};

//==============================================================================
