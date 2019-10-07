#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>

class Mesh
{
	public:

	int n;
	std::vector<double> quadSet;
	std::vector<double> mu;
	std::vector< std::vector<double> > ordinates;
        Mesh();
	void calcMu();
	void calcQuadSet();
	void printQuadSet();
};
#endif
