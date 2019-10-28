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
        StartingAngle(Mesh::Mesh myMesh,YAML::Node myInput);

        private:
        // private functions
        YAML::Node * input;
        Mesh::Mesh * mesh;
};

//==============================================================================
