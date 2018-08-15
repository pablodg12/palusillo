#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "palu.h"


/*Estas funciones se explican en palu.h*/

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
	for(int h = 0; h < size; h++){
		aux[h] = A[i*size + h];
	}
	for(int h = 0; h < size; h++){
		A[i*size + h] = A[(j-i) + h];
	}
	for(int h = 0; h < size; h++){
		A[(j-i) + h] = aux[h];
	}
	return(A);
}

double* solve_ALU(int size, double* A, double* b){
	double* c =  malloc(sizeof(double)*size);
	c[0] = b[0];
	for(int i = 1;i<size;i++){
		float tmp = 0;
		for(int z=0;z<i;z++){
			tmp = tmp + A[i*size + z] * c[z];
		}
		c[i] = (1.0)* (b[i]-tmp);
	}
	b[size-1] = (1.0/A[size*size-1]) * c[size-1];
	for(int i = size-2;i>-1;i--){
		float tmp = 0;
		for(int z=i+1;z<size;z++){
			tmp = tmp + A[i*size + z] * b[z];
		}
		b[i] = (1.0/A[i*size+i])* (c[i]-tmp);
		}
		free(c);
		return(b);
	}


/* ############################################################*/

double* solve_palu(double* A, double* b,int size){
	int p_index = 0;
	double tmp = 0;

	/*Ciclo de PALU, se realizan las permutaciones, junto a la reducción para formar 
	la matriz superior e inferior */

	for(int j=0;j<size-1;j++){
		p_index = argmax(size,A,j);
		if(p_index > 0){
			A = row_perm(size,A,j,p_index);
			p_index = (int) p_index/size;
			tmp = b[j];
			b[j] = b[p_index];
			b[p_index]=tmp;
		}

		for(int i = j+1; i < size; i++){
			A[i*size + j] = A[i*size + j]/A[j*size + j];
			for(int h = j+1; h < size; h++){
				A[i*size+h] = A[i*size+h] - A[i*size + j]*A[j*size + h];
			}
		}

	}

	b = solve_ALU(size,A,b);

	return(b);
};
