#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>

using namespace std;

class Mesh
{
	public:
	Mesh();  	
  	int n;		
  	vector< vector<double> > quadSet;
  	vector< vector<double> > alpha;
  	void calcQuadSet();
	void printQuadSet();
	
        private:
	vector<double> mu;
  	vector< vector<double> > ordinates;
  	void calcMu();	
	void calcAlpha();
	int quad_index(int p,int q);
	int low_quad_index(int p,int q);
};
#endif
