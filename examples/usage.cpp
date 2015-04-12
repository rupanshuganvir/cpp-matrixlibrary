#include <iostream>
#include "matrix.h"

using namespace std;

int main(){

	double A[]= {2,1,1,4,-6,0,-2,7,2};
	matrix<double> P(3,3,A);
	cout<<P<<endl;

	cout<<reduced(P)<<endl;

	cout<<C_space(P)<<endl;

	return 0;
}
