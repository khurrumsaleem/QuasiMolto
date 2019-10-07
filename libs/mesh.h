#ifndef MESH_H
#define MESH_H

#include <iostream>

class Mesh
{
	public:

	int n;
	double* quadSet;
	double* mu;
	double* ordinates;
        Mesh();
	void calcMu();
	void calcQuadSet();
	void printQuadSet();
};
#endif
