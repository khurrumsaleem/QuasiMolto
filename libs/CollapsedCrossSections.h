#ifndef COLLAPSEDCROSSSECTIONS_H
#define COLLAPSEDCROSSSECTIONS_H

#include "Mesh.h"

using namespace std; 
using namespace arma;

//==============================================================================
//! Class the holds collapsed one-group nuclear data

class CollapsedCrossSections
{
        public:
        CollapsedCrossSections(int nZ,int nR);
	Eigen::MatrixXd sigT,sigS,sigF,chiP,chiD,neutV,nu; 
};

//==============================================================================

#endif
