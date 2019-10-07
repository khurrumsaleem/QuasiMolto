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

	void calcMu();
	void calcQuadSet();
	void printQuadSet();
};
#endif
