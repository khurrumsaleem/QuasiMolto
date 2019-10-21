// test.cpp

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void printTransport() {

     field<vec> flux;
     flux.set_size(10,5);
     
     for (int irows = 0; irows < flux.n_rows; ++irows) {
     	for (int icols = 0; icols < flux.n_cols; ++icols) {
		flux(irows,icols) = randu<vec>(2,1);
     	}
     }

     cout << "Printed from transport" << endl;	

}
