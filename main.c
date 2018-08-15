#include <stdio.h>
#include <stdlib.h>	
#include "palu.h"


int main( int argc, const char* argv[] )
{
	int size = 6;
	double A[ ] = {3.0,-1.0,0.0,0.0,0.0,0.5,-1.0,3.0,-1.0,0.0,0.5,0.0,0.0,-1.0,3.0,-1.0,0.0,0.0,0.0,0.0,-1.0,3.0,-1.0,0.0,0.0,0.5,0.0,-1.0,3.0,-1.0,0.5,0.0,0.0,0.0,-1.0,3.0};
	double b[ ] = {2.5,1.5,1.0,1.0,1.5,2.5};
	double* x;
	x = palu_decomp(A,b,size);
	for (int i = 0; i < size; ++i)
	{
		printf("%f\n", x[i]);
	}
	return 0;
}
