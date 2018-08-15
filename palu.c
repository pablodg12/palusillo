#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "palu.h"


/*Estas funciones se explican en palu.h*/

double* generate_identity(int size){
	double* aux =  malloc(sizeof(double)*size*size);
	memset(aux,0,sizeof(double)*size*size);
	for(int i=0;i<size;i++){
		aux[i + i*size] = 1;
	}
	return(aux);
};

double* generate_zeros(int size){
	double* aux =  malloc(sizeof(double)*size*size);
	memset(aux,0,sizeof(double)*size*size);
	return(aux);
};

int argmax(int size, double* A,int j){
	int max_index = -1;
	float max_value = -1.0;
	for(int i=0;i<size;i++){
		if(i>=j){
			if(fabs(A[(j)+size*i])>max_value){
					max_value = fabs(A[(j)+size*i]);
					max_index = (j)+size*i;
				}
			}
		}
	return(max_index);
};

double* row_perm(int size, double* A, int i, int j){
	double aux[size];
	double* U = A;
	for(int h = 0; h < size; h++){
		aux[h] = U[i*size + h];
	}
	for(int h = 0; h < size; h++){
		U[i*size + h] = U[(j-i) + h];
	}
	for(int h = 0; h < size; h++){
		U[(j-i) + h] = aux[h];
	}
	return(U);
}

double* dot(int size, double* P, double* b){
	double* aux =  malloc(sizeof(double)*size);
	for(int i=0;i<size;i++){
		float tmp = 0;
		for(int j=0;j<size;j++){
			tmp = tmp + P[i*size + j] *b[j];
		}
		aux[i] = tmp;
	}
	return(aux);
}

double* solve_triangular(int size, int upper, double* A, double* b){
	double* x =  malloc(sizeof(double)*size);
	if(upper == 1){
		x[size-1] = (1.0/A[size*size-1]) * b[size-1];
		for(int i = size-2;i>-1;i--){
			float tmp = 0;
			for(int z=i+1;z<size;z++){
				tmp = tmp + A[i*size + z] * x[z];
			}
			x[i] = (1.0/A[i*size+i])* (b[i]-tmp);
		}
		return(x);
	}
	else{
		x[0] = (1.0/A[0]) * b[0];
		for(int i = 1;i<size;i++){
			float tmp = 0;
			for(int z=0;z<i;z++){
				tmp = tmp + A[i*size + z] * x[z];
			}
			x[i] = (1.0/A[i*size+i])* (b[i]-tmp);
		}
		return(x);
	}
}

/* ############################################################*/

double* palu_decomp(double* A, double* b,int size){
	int p_index = 0;
	double* P = generate_identity(size);
	double* L = generate_zeros(size);
	double* U = A;
	double* c;
	double* x;

	/*Ciclo de PALU, se realizan las permutaciones, junto a la reducción para formar 
	la matriz superior e inferior */

	for(int j=0;j<size-1;j++){
		p_index = argmax(size,U,j);
		if(p_index > 0){
			P = row_perm(size,P,j,p_index);
			L = row_perm(size,L,j,p_index);
			U = row_perm(size,U,j,p_index);
		}

		for(int i = j+1; i < size; i++){
			L[i*size + j] = U[i*size + j]/U[j*size + j];
			for(int h = 0; h < size; h++){
				U[i*size+h] = U[i*size+h] - L[i*size + j]*U[j*size + h];
			}
		}

	}
	for(int i=0;i<size;i++){
		L[i*size+i] = 1;
	}

	/*Print de las matrices P,L,U*/
	printf("---------P--------\n");
	for(int ii = 0; ii < size; ii++){
		for(int jj = 0; jj < size; jj++){
			printf("%f ",P[ii*size + jj]);
		}
		printf("\n");
	}
	printf("---------L---------\n");
	for(int ii = 0; ii < size; ii++){
		for(int jj = 0; jj < size; jj++){
			printf("%f ",L[ii*size + jj]);
		}
		printf("\n");
	}
	printf("---------U-------\n");
	for(int ii = 0; ii < size; ii++){
		for(int jj = 0; jj < size; jj++){
			printf("%f ",U[ii*size + jj]);
		}
		printf("\n");
	}
	printf("----------------\n");

	/* Solve dot(P,b) junto a la triangular superior e inferior*/
	b = dot(size, P, b);
	printf("---------b-------\n");
	for (int i = 0; i < size; ++i)
	{
		printf("%f\n", b[i]);
	}
	c = solve_triangular(size,0,L,b);
	printf("---------c-------\n");
	for (int i = 0; i < size; ++i)
	{
		printf("%f\n", c[i]);
	}
	x = solve_triangular(size,1,U,c);
	printf("---------x-------\n");

	/*Se libera la memoría pedida */
	free(P);
	free(L);
	free(b);
	free(c);
	return(x);

};
