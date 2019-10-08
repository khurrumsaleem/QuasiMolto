#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>

class Mesh
{
	public:

	int n;
	std::vector< std::vector<double> > quadSet;
	std::vector<double> mu;
	std::vector< std::vector<double> > ordinates;
	std::vector< std::vector<double> > alpha;
        Mesh();
	void calcMu();
	void calcQuadSet();
	void calcAlpha();
	void printQuadSet();
	int quad_index(int p, int q);
	int low_quad_index(int p,int q);
};
#endif
